#include "dydx_function.h"

sp_cx_mat dydxFunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy)
{

	int ng = nx * ny;

	cx_vec a = p % join_cols(cx_vec(nx + 1), r.rows(0, ng - (nx + 1) - 1)) % join_cols(cx_vec(nx), 1 / q.rows(0, ng - nx - 1));
	cx_vec b = -p % join_cols(cx_vec(nx - 1), r.rows(0, ng - (nx - 1) - 1)) % join_cols(cx_vec(nx), 1 / q.rows(0, ng - nx - 1));
	cx_vec c = -p % join_cols(r.rows(nx - 1, ng - 1), cx_vec(nx - 1)) % join_cols(1 / q.rows(nx, ng - 1), cx_vec(nx));
	cx_vec d = p % join_cols(r.rows(nx + 1, ng - 1), cx_vec(nx + 1)) % join_cols(1 / q.rows(nx, ng - 1), cx_vec(nx));



	for (int i = 0;i < ng;i += nx) {
		a(i) = 0;
		c(i) = 0;

	}
	for (int i = nx - 1;i < ng;i += nx) {
		b(i) = 0;
		d(i) = 0;
	}

	sp_cx_mat s = spdiags(join_rows(d, c, b, a) / 4 / dx / dy,
		ivec{ -nx - 1, -nx + 1, nx - 1, nx + 1 }, ng, ng).st();

	return s;
}

