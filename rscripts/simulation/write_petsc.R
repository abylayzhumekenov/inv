# ------------------------------------------------------------------------------
# ---------- FUNCTIONS FOR BINARY WRITE ----------------------------------------
# ------------------------------------------------------------------------------

# matrix
write_petsc_mat = function(Q, filename){
    # encode the matrix
    Q = as(Q, "RsparseMatrix")
    x = list(classid = 1211216,
             nrows = nrow(Q),
             ncols = ncol(Q),
             nnz = length(Q@x),
             nnz_row = diff(Q@p),
             nnz_i = Q@j,
             nnz_val = Q@x)
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$ncols), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_row), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nnz_i), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$nnz_val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}

# vector
write_petsc_vec = function(y, filename){
    # encode the vector
    x = list(classid = 1211214,
             nrows = length(y),
             val = drop(y))
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.double(x$val), con = fwrite, size = 8, endian = "swap")
    close(con = fwrite)
}

# index set
write_petsc_is = function(is, filename){
    # encode the index set
    x = list(classid = 1211218,
             nrows = length(is),
             val = drop(is)-1)
    
    # write the binary data
    fwrite = file(filename, "wb")
    writeBin(object = as.integer(x$classid), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$nrows), con = fwrite, size = 4, endian = "swap")
    writeBin(object = as.integer(x$val), con = fwrite, size = 4, endian = "swap")
    close(con = fwrite)
}
