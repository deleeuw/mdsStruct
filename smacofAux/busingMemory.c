double **getmatrix(const size_t nr, const size_t nc, const double c)
// allocates matrix space on the heap
{
  double **ptr = 0;
  if (nr == 0 || nc == 0)
    return ptr;
  double *block = 0;
  ptr = (double **)calloc(nr, sizeof(double *))
  block = (double *)calloc(nr * nc, sizeof(double));
  ptr--;
  block--;
  for (size_t i = 1, im1 = 0; i <= nr; i++, im1++) {
    ptr[i] = &block[im1 * nc];
    for (size_t j = 1; j <= nc; j++)
      ptr[i][j] = c;
  }
  return ptr;
} // getmatrix

void freematrix(double **a)
// de-allocates matrix space from the heap
{
  free(++a[1]);
  free(++a);
} // freematrix