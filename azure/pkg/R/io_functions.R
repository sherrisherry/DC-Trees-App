keycache <- read.csv('~/vars/azcodes.csv', header = T, stringsAsFactors = F) # the database of credentials

nv_cols <- function(nm, cl, obj, io = TRUE){
  switch(obj,
         bridge = {# "whem","mena","ssa","deur","asia" in imf_reg
           cols <- c(rep('integer',7), rep('character',9))
           names(cols) <- c('un_code','un_start','un_end','wb_code','imf_code','d_gfi','d_dev',
                            'imf_reg','imf_nm','un_nm','un_nm_en_full','un_nm_en_abbr','un_note','iso2','iso3','wb_nm')
         },
         geo = {
           cols <- c(rep("integer",4),"numeric",rep("integer",4))
           names(cols) <- c("j","i","area_j","area_i","distw","d_landlocked_j","d_landlocked_i","d_contig","d_conti")
         },
         eia = {
           cols <- rep("integer",4)
           names(cols) <- c("t","i","j","d_rta")
         },
         hkrx = {
           cols <- c(rep("integer",3),"character","numeric",rep("integer",2),"numeric")
           names(cols) <- c("t","origin_hk","consig_hk","k","vrx_hkd","origin_un","consig_un","vrx_un")
         })
  if(missing(nm))nm <- names(cols)
  if(!missing(cl))cols[nm] <- cl
  if(io)cols[setdiff(names(cols), nm)] <- 'NULL' else cols[nm] <- 'NULL'
  return(cols)
}

in_bridge <- function(nm, cl, logf, max_try = 10, io = TRUE){
  cols <- nv_cols(nm, cl, 'bridge', io)
  pth <- deco_path_stor('dlake02/coding-systems/bridge.csv')
  batchscr::ecycle(bridge <- az_read_blob(FUN = function(x)read.csv(x, colClasses=cols, header=T, fileEncoding = 'UTF8'),
                                          object = pth$obj, container = az_container(pth)),
                   {if(!missing(logf))logf(paste('0000', '!', 'loading bridge.csv failed', sep = '\t')); return(NULL)}, max_try)
  return(bridge)
}

in_geo <- function(nm, cl, logf, max_try = 10, io = TRUE){
  cols <- nv_cols(nm, cl, 'geo', io)
  pth <- deco_path_stor('cachetemp/supplemental/CEPII_GeoDist.csv')
  batchscr::ecycle(geo <- az_read_blob(FUN = function(x)read.csv(x, colClasses=cols, header=T, na.strings=''),
                                       object = pth$obj, container = az_container(pth)),
                   {if(!missing(logf))logf(paste('0000', '!', 'loading CEPII_GeoDist.csv failed', sep = '\t')); return(NULL)}, max_try)
  return(geo)
}

in_eia <- function(nm, cl, logf, max_try = 10, io = TRUE){
  cols <- nv_cols(nm, cl, 'eia', io)
  pth <- deco_path_stor('cachetemp/supplemental/EIA.csv.bz2')
  tmp <- tempfile()
  batchscr::ecycle(AzureStor::storage_download(src = pth$obj, container = az_container(pth), dest = tmp, overwrite = TRUE),
                   {if(!missing(logf))logf(paste('0000', '!', 'retrieving EIA file failed', sep = '\t')); return(NULL)}, max_try)
  batchscr::ecycle(eia <- read.csv(bzfile(tmp), header=T, colClasses=cols, na.strings="", stringsAsFactors = F),
                   {if(!missing(logf))logf(paste('0000', '!', 'loading file failed', sep = '\t')); return(NULL)},
                   max_try, cond = is.data.frame(eia) && nrow(eia)>10)
  return(eia)
}

in_hkrx <- function(yr, nm, cl, logf, max_try = 10, io = TRUE){
  cols <- nv_cols(nm, cl, 'hkrx', io)
  pth <- deco_path_stor('cachetemp/supplemental')
  tmp <- tempfile()
  batchscr::ecycle(AzureStor::storage_download(src = paste('HK', yr, 'rx.csv.bz2', sep = '_'), container = az_container(pth), dest = tmp, overwrite = TRUE),
                   {if(!missing(logf))logf(paste(yr, '!', 'retrieving hkrx file failed', sep = '\t')); return(NULL)}, max_try)
  batchscr::ecycle(hk <- read.csv(bzfile(tmp), header=T, colClasses=cols, na.strings="", stringsAsFactors = F),
                   {if(!missing(logf))logf(paste(yr, '!', 'loading hkrx file failed', sep = '\t')); return(NULL)},
                   max_try, cond = is.data.frame(hk) && nrow(hk)>10)
  return(hk)
}

# opt can be missing, dev, adv, or a region abbrivation
gfi_cty <- function(opt, logf, max_try = 10){
  cols <- if(!missing(opt))switch(opt, dev = 'd_dev', adv = c('d_gfi', 'd_dev'), c('d_gfi', 'imf_reg'))else 'd_gfi'
  bridge <- in_bridge(c(cols, 'un_code'), logf = logf, max_try = max_try)
  if(is.null(bridge))return(NULL)
  bridge <- unique(bridge)
  cty <- if(length(cols)==1)subset(bridge, bridge[, cols]==1, 'un_code', drop = T)
  else switch(class(bridge[, setdiff(cols, 'd_gfi')]),
              integer = subset(bridge, abs(bridge[, cols[1]]-bridge[, cols[2]])==1, 'un_code', drop = T),
              character = bridge$un_code[bridge$d_gfi==1 & bridge$imf_reg==opt])
  if(!missing(logf))logf(paste('0000', ':', 'decided cty', sep = '\t'))
  return(cty)
}

deco_path_stor <- function(path){
  vct <- strsplit(path, '/')[[1]]
  path <- list()
  path$acct <- vct[1]; path$cont <- vct[2]
  path$obj <- ifelse(length(vct)>2, paste(vct[-1:-2], sep = '/'), '')
  return(path)
}

az_ep <- function(resrc, serv){
  switch(serv,
         blob = {prefix = "https://"; suffix = ".blob.core.windows.net/"},
         dlake = {prefix = "https://"; suffix = ".dfs.core.windows.net/"},
         file = {prefix = "https://"; suffix = ".file.core.windows.net/"},
         queue = {prefix = "https://"; suffix = ".queue.core.windows.net/"})
  return(paste(prefix, resrc, suffix, sep = ""))
}

az_key <- function(resrc, serv){
  return(keycache$key[keycache$name==resrc & keycache$service==serv])
}

az_container <- function(obj){
  ep <- AzureStor::storage_endpoint(az_ep(obj$acct, 'blob'), key = az_key(obj$acct, 'storage'))
  return(AzureStor::storage_container(ep, obj$cont))
}

az_read_blob <- function(FUN, object, container){
  tmp <- tempfile()
  AzureStor::storage_download(container, object, tmp)
  blob <- FUN(tmp)
  unlink(tmp)
  return(blob)
}

az_write_blob <- function(obj, FUN, object, container){
  tmp <- tempfile()
  FUN(obj, tmp)
  AzureStor::storage_upload(container, tmp, object)
}
