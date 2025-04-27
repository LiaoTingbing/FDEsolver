#include "dxdx_function.h"

sp_cx_mat dxdxFunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy)
{
	int ng = nx * ny;
	cx_vec a = p % join_cols(cx_vec(1), r.rows(0, ng - 2))
		/ (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q);

	cx_vec b = -p % r % (
		1 / (join_cols(cx_vec(1), q.rows(0, ng - 2)) + q)
		+ 1 / (q + join_cols(q.rows(1, ng - 1), cx_vec(1)))
		);
	cx_vec c = p % join_cols(r.rows(1, ng - 1), cx_vec(1))
		/ (q + join_cols(q.rows(1, ng - 1), cx_vec(1)));

	for (int i = 0;i < ng;i += nx)
		a(i) = 0;

	for (int j = nx - 1;j < ng;j += nx)
		c(j) = 0;

	sp_cx_mat s = spdiags(join_rows(c, b, a) * 2 / dx / dx,
		ivec{ -1,0,1 }, ng, ng).st();
	return s;

}

sp_mat dxdxFuncReal(const vec& p, vec& q, vec& r, int nx, int ny, double dx, double dy)
{
	int ng = nx * ny;
	vec a = p % join_cols(vec(1), r.rows(0, ng - 2))
		/ (join_cols(vec(1), q.rows(0, ng - 2)) + q);

	vec b = -p % r % (
		1 / (join_cols(vec(1), q.rows(0, ng - 2)) + q)
		+ 1 / (q + join_cols(q.rows(1, ng - 1), vec(1)))
		);
	vec c = p % join_cols(r.rows(1, ng - 1), vec(1))
		/ (q + join_cols(q.rows(1, ng - 1), vec(1)));

	for (int i = 0;i < ng;i += nx)
		a(i) = 0;

	for (int j = nx - 1;j < ng;j += nx)
		c(j) = 0;
	sp_mat s = spdiags(join_rows(c, b, a) * 2 / dx / dx,
		ivec{ -1,0,1 }, ng, ng).st();
	return s;
}
