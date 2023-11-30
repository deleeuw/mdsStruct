/* HS, p 230 */

int main(void) {
	int n= 5, m = 3;
	double rect[n][m]; /* VLA cannot be initialized */
	double (*p)[m]; /* p is an object of type double[m] */
	p = rect; /* same as p = &rect[0], pointer to first row of rect */
	p++; /* now p = &rect[1], second row of rect */
}