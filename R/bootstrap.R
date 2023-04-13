## bootstraping
#function for blockbootstapping. Returns the reshuffled row ids.
#nrows_df: nomber of rows of the original precipitation matrix
#block_size: size of block, eg 7 for 1 week in the case of daily data
# returns new "row id" for a block bootstap.
block_boot_function = function(nrows_df, block_size){

  ndays = nrows_df#nrow(new)
  block_size = block_size
  nblocks = ceiling(ndays/block_size)

  row_id = 1:(nblocks*block_size)


  # convert to matrix
  row_id_mat = matrix(row_id, nrow = nblocks, ncol = block_size, byrow = T)

  #sample with replacment , the nomber of blcoks
  bootstap_blocks = sample(1:nblocks, size = nblocks, replace = T)
  #extract the bootrap matrxi
  bootstap_matrix = row_id_mat[bootstap_blocks,]
  #convert to vector
  bootstap_row_id = as.vector(t(bootstap_matrix))[1:ndays] #1:ndays to enusure the same
  #length as the original matrix, that means the last block might not be of size "block size' if ndays/block_size
  #is not a whole number
  return(bootstap_row_id)
}
