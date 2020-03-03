library(tidyverse)
library(SparseM)
library(gridExtra)
library(base2grob)

#----

bits = 5

plate_df = expand.grid(c('a','b','c','d','e','f','g','h'), 1:12)
# plate_arr = array(0, dim=c(8,16,bits))

split_bit_into_n_blocks = function(code_number, n_blocks) {
  
  run_length = code_number/n_blocks/2 # divide by 2 to account for label AND unlabelled wells
  bit_value = 0
  code_vec = NA
  for (i in 1:code_number) {
    # if the current position is divisible by the run_length, switch bit value
    if ((i-1) %% run_length == 0) {
      bit_value = ifelse(bit_value == 1, 0, 1)
    }
    code_vec[i] = bit_value
  }
  return(code_vec)
}

split_bit_into_n_blocks(32, 2)

rle(split_bit_into_n_blocks(32, 4))

construct_code_matrix = function(bit_count) {
  
  code_number = 2**bit_count
  mat = matrix(0, nrow=bit_count, ncol=code_number)
  
  for (bit_position in 1:nrow(mat)) {
    
    if (bit_position == 1) {
      block_number = 1
    } else {
      block_number = block_number * 2
    }
    
    code_vec = split_bit_into_n_blocks(code_number, block_number)
    mat[bit_position,] = code_vec
    
  }
  
  return(mat)
  
}

code_mat = construct_code_matrix(bit_count=5)

construct_code_array = function(bit_count, plate_size=96, plate_cols=12) {
  
  plate_rows = plate_size/plate_cols
  code_number = 2**bit_count
  
  code_mat = construct_code_matrix(bit_count) 
  
  code_cols = code_number/plate_rows
  code_arr = array(t(code_mat), dim=c(plate_rows, code_cols, bit_count))
  
  return(code_arr)
  
}

code_arr = construct_code_array(bit_count=5)
code_arr

SparseM::image(as.matrix.csr(code_arr[,,5]), main="Dye 5")

plot_list = list()
for (dye_index in 1:dim(code_arr)[3]) {
  plot_list[[dye_index]] = base2grob(~SparseM::image(as.matrix.csr(code_arr[,,dye_index]), main=paste("Dye", dye_index, sep=": ")))
}

plot_grid = grid.arrange(
  grobs = plot_list,
  nrow=3
)

ggsave(plot_grid, filename="five-bit_dye_positions.pdf")

paste_code = function(code_vec) {
  return(paste0(code_vec, collapse=''))
}

assemble_codebook = function(code_array) {
  return(apply(code_array, paste_code, MARGIN=1:2))
}

codebook = assemble_codebook(code_arr)
codebook_tib = as.tibble(codebook) %>%
  mutate(row=toupper(letters[1:8]))
names(codebook_tib)[1:4] = 1:4
codebook_tib = codebook_tib %>%
  gather(key=column, val=code, -row) %>%
  mutate(column = as.numeric(column))

plate_layout = read_tsv('20200303_strain_plate_positions.tsv') %>%
  mutate(library = ifelse(`Plate Column` >= 5, 2, 1),
         column = ifelse(library == 1, `Plate Column`, `Plate Column` - 4)) %>%
  dplyr::rename(row = `Plate Row`) %>%
  left_join(codebook_tib, by=c("row"="row", "column"="column")) %>%
  dplyr::rename(BKK_tag=`Duplicate minus misc_RNA/missing library data`,
                gene=`gene deletion`) %>%
  dplyr::select(BKK_tag, gene, row, column, code, library) %>%
  arrange(row,column)
plate_layout

plate_layout %>%
  write.csv(file='20200303_vanco_and_cipro_codebook.csv', row.names=FALSE)



code_arr = construct_code_array(bit_count=5)
plot_list = list()
for (dye_index in 1:dim(code_arr)[3]) {
  plot_list[[dye_index]] = base2grob(~SparseM::image(as.matrix.csr(code_arr[,,dye_index]), main=paste("Dye", dye_index, sep=": ")))
}

plot_grid = grid.arrange(
  grobs = plot_list,
  nrow=3
)

ggsave(plot_grid, filename="five-bit_dye_positions.pdf")

#----
make_code_strs = function(x) {
  return(paste(x, collapse=''))
}

codes = apply(code_arr, FUN=make_code_strs, MARGIN=1:2)

any(duplicated(as.vector(codes)))

codes


