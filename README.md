# NRSTDP

### Program archive for article 'Neurodynamical Role of STDP in Storage and Retrieval of Associative Information'.

1.  File 'image_group.m' is for the simulation of 'retrieval of grouped images' section in the article.
    Basically performs association of five 32 x 32 pixels of orchestral instrument images used in the article.
    To raw run this file, please include the following 6 images files 'violin.png', 'trumpet.png', 'harp.png', 'piano.png', 'timpani.png', and 'forest.png', and 'gsprocess.m' in the folder containing this file.
    User can arbitrarily replace given images and parameters such as the dimension of the images, etc.

2.  File 'composite_struct.m' is for the simulation of 'multiple groups of memory with composite structure' section in the article.
    Basically uses the same sentences S1, S2, and S3 in the article.
    Please include 'gsprocess.m' in the folder containing this file.
    User can randomly replace the structure of given sentences and perform the same tasks of retrieval.

3.  File 'gsprocess.m' is a simple function file conducting Gram-Schmidt process for the set of column vectors of an arbitrary matrix.
