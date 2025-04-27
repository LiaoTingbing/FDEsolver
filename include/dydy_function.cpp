#include "dydy_function.h"

sp_cx_mat dydyFunc(const cx_vec& p, const cx_vec& q, const cx_vec& r, int nx, int ny, double dx, double dy)
{
	int ng = nx * ny;

	cx_vec a = p % join_cols(cx_vec(nx), r.rows(0, ng - nx - 1))
		/ (join_cols(cx_vec(nx), q.rows(0, ng - nx - 1)) + q);
	cx_vec b = -p % r % (1 / (join_cols(cx_vec(nx), q.rows(0, ng - nx - 1)) + q)
		+ 1 / (q + join_cols(q.rows(nx, ng - 1), cx_vec(nx))));
	cx_vec c = p % join_cols(r.rows(nx, ng - 1), cx_vec(nx))
		/ (q + join_cols(q.rows(nx, ng - 1), cx_vec(nx)));

	sp_cx_mat s = spdiags(join_rows(c, b, a) * 2 / dy / dy,
		ivec{ -nx,0,nx }, ng, ng).st();
	return s;



}

sp_mat dydyFuncReal(const vec& p, const vec& q, const vec& r, int nx, int ny, double dx, double dy)
{
	int ng = nx * ny;

	vec a = p % join_cols(vec(nx), r.rows(0, ng - nx - 1))
		/ (join_cols(vec(nx), q.rows(0, ng - nx - 1)) + q);
	vec b = -p % r % (1 / (join_cols(vec(nx), q.rows(0, ng - nx - 1)) + q)
		+ 1 / (q + join_cols(q.rows(nx, ng - 1), vec(nx))));
	vec c = p % join_cols(r.rows(nx, ng - 1), vec(nx))
		/ (q + join_cols(q.rows(nx, ng - 1), vec(nx)));

	sp_mat s = spdiags(join_rows(c, b, a) * 2 / dy / dy,
		ivec{ -nx,0,nx }, ng, ng).st();
	return s;


}
