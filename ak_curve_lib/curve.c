#include "curve.h"
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
void ak_point_set_point(ak_point ep1, ak_point ep2, ak_curve ec)
{
    memcpy(ep1->x, ep2->x, ec->size * sizeof(ak_uint64));
    memcpy(ep1->y, ep2->y, ec->size * sizeof(ak_uint64));
    memcpy(ep1->z, ep2->z, ec->size * sizeof(ak_uint64));
}

void ak_point_set(ak_point ep, ak_curve ec) {
    memcpy(ep->x, ec->point.x, ec->size * sizeof(ak_uint64));
    memcpy(ep->y, ec->point.y, ec->size * sizeof(ak_uint64));
    memcpy(ep->z, ec->point.z, ec->size * sizeof(ak_uint64));
}

void ak_epoint_set_as_unit(ak_point ep, ak_curve ec)
{
    ak_mpzn_set_ui(ep->x, ec->size, 0);
    ak_mpzn_set_ui(ep->y, ec->size, 1);
    ak_mpzn_set_ui(ep->z, ec->size, 1);
}
void ak_point_set_as_unit(ak_point ep, ak_curve ec){
    if (ec->cform==0) ak_wpoint_set_as_unit((ak_wpoint)ep,(ak_wcurve)ec);
    if (ec->cform==1)   ak_epoint_set_as_unit(ep,ec);
}
void ak_point_reduce(ak_point ep, ak_curve ec)
{
    if(ec->cform==0) ak_wpoint_reduce((ak_wpoint)ep,(ak_wcurve)ec);
    if(ec->cform==1){
        if (ak_mpzn_cmp_ui(ep->z, ec->size, 0)) return;
	ak_mpznmax u, one = ak_mpznmax_one;
	ak_mpzn_set_ui(u, ec->size, 2);
	ak_mpzn_sub(u, ec->p, u, ec->size);
    ak_mpzn_modpow_montgomery(u, ep->z, u, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(u, u, one, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->x, ep->x, u, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->y, ep->y, u, ec->p, ec->n, ec->size);
    ak_mpzn_set_ui(ep->z, ec->size, 1);
    }
}

bool_t ak_epoint_is_ok(ak_point ep, ak_curve ec) {
	ak_mpznmax a, b, c, d;
    ak_mpzn_mul_montgomery(a, ep->x, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(d, a, ec->a, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, ep->y, ep->y, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(d, d, b, ec->p, ec->size);
    ak_mpzn_mul_montgomery(c, ep->z, ep->z, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(d, d, c, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(c, c, c, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(a, a, b, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(a, a, ec->b, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(a, a, c, ec->p, ec->size);
	if (ak_mpzn_cmp(a, d, ec->size)) return ak_false;
	return ak_true;
}
bool_t ak_point_is_ok(ak_point ep, ak_curve ec) {
    if (ec->cform==0) return ak_wpoint_is_ok((ak_wpoint)ep,(ak_wcurve)ec);
    if (ec->cform==1) return ak_epoint_is_ok(ep,ec);
    return ak_false;
}
void ak_epoint_double(ak_point ep, ak_curve ec) {
	ak_mpznmax b, c, d, e, f, h, j, st;
    ak_mpzn_add_montgomery(b, ep->x, ep->y, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, b, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(c, ep->x, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(d, ep->y, ep->y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(e, ec->a, c, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(f, e, d, ec->p, ec->size);
    ak_mpzn_mul_montgomery(h, ep->z, ep->z, ec->p, ec->n, ec->size);
	ak_mpzn_lshift_montgomery(j, h, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, j, ec->size);
	ak_mpzn_add_montgomery(j, f, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, c, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->x, b, j, ec->p, ec->n, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(e, e, st, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->y, f, e, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->z, f, j, ec->p, ec->n, ec->size);
}
void ak_point_double(ak_point ep, ak_curve ec){
    if (ec->cform==0) ak_wpoint_double((ak_wpoint)ep,(ak_wcurve)ec);
    if (ec->cform==1)   ak_epoint_double(ep,ec);
}
void ak_epoint_add(ak_point ep1, ak_point ep2, ak_curve ec) {
	ak_mpznmax a, b, c, d, e, f, g, st;
    ak_mpzn_mul_montgomery(a, ep1->z, ep2->z, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(b, a, a, ec->p, ec->n, ec->size); 
    ak_mpzn_mul_montgomery(c, ep1->x, ep2->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(d, ep1->y, ep2->y, ec->p, ec->n, ec->size);
	ak_mpzn_mul_montgomery(e, c, d, ec->p, ec->n, ec->size); 
    ak_mpzn_mul_montgomery(e, e, ec->b, ec->p, ec->n, ec->size);
	ak_mpzn_sub(st, ec->p, e, ec->size);
	ak_mpzn_add_montgomery(f, b, st, ec->p, ec->size);
	ak_mpzn_add_montgomery(g, b, e, ec->p, ec->size); 
    ak_mpzn_add_montgomery(e, ep1->x, ep1->y, ec->p, ec->size);
    ak_mpzn_add_montgomery(b, ep2->x, ep2->y, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, e, ec->p, ec->n, ec->size); 
	ak_mpzn_sub(st, ec->p, c, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_sub(st, ec->p, d, ec->size);
	ak_mpzn_add_montgomery(b, b, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, f, ec->p, ec->n, ec->size); 
    ak_mpzn_mul_montgomery(ep1->x, b, a, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, ec->a, c, ec->p, ec->n, ec->size);
	ak_mpzn_sub(st, ec->p, b, ec->size);
	ak_mpzn_add_montgomery(b, d, st, ec->p, ec->size);
	ak_mpzn_mul_montgomery(b, b, g, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep1->y, b, a, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep1->z, f, g, ec->p, ec->n, ec->size);
}
void ak_point_add(ak_point ep1, ak_point ep2, ak_curve ec){
    if (ec->cform==0) ak_wpoint_add((ak_wpoint)ep1,(ak_wpoint)ep2,(ak_wcurve)ec);
    if (ec->cform==1)   ak_epoint_add(ep1,ep2,ec);
}

void ak_point_pow(ak_point ep1, ak_point ep2, ak_uint64 * a, size_t size, ak_curve ec) {
    ak_uint64 uk;
	long long int i, j;
    struct point Q, R;
    ak_point_set_as_unit(&Q, ec);
    ak_point_set_point(&R, ep2, ec);
	for (i = size - 1; i >= 0; i--) {
		uk = a[i];
		for (j = 0; j < 64; j++) {
            if (uk & 0x8000000000000000LL) {ak_point_add(&Q, &R, ec); ak_point_double(&R, ec);}
            else {	ak_point_add(&R, &Q, ec); ak_point_double(&Q, ec);}
			uk <<= 1;
		}
	}
    ak_point_set_point(ep1, &Q, ec);
}

void ak_set_curve_form_e(ak_curve ec){
    ak_mpznmax a,b;
    memcpy(a, ec->edw.a, ec->size * sizeof(ak_uint64));
    memcpy(b, ec->edw.b, ec->size * sizeof(ak_uint64));
    memcpy(ec->a, ec->edw.a, ec->size * sizeof(ak_uint64));
    memcpy(ec->b, ec->edw.b, ec->size * sizeof(ak_uint64));
    memcpy(ec->a, a, ec->size * sizeof(ak_uint64));
    memcpy(ec->b, b, ec->size * sizeof(ak_uint64));
    ec->cform=ec->cform==1?0:1;
}
void ak_epoint_to_wpoint(ak_point wp,ak_point ep, ak_curve ec) {
    if (ak_mpzn_cmp_ui(ep->z, ec->size, 0)) {
        ak_set_curve_form_e(ec);
        ak_point_set_as_unit(wp, ec);
		return;
	}
	ak_mpznmax a_;
    ak_mpzn_sub(wp->z, ec->p, ep->y, ec->size);
    ak_mpzn_add_montgomery(wp->z, wp->z, ep->z, ec->p, ec->size);
    ak_mpzn_mul_montgomery(wp->z, wp->z, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(a_, ep->y, ep->z, ec->p, ec->size);
    ak_mpzn_mul_montgomery(a_, a_, ec->edw.s, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(wp->y, a_, ep->z, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(wp->x, a_, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(a_, ec->edw.t, wp->z, ec->p, ec->n, ec->size);
	ak_mpzn_add_montgomery(wp->x,wp->x, a_, ec->p, ec->size);
    ak_set_curve_form_e(ec);
}
void ak_wpoint_to_epoint(ak_point ep,ak_point wp,ak_curve ec) {
    if (ak_mpzn_cmp_ui(wp->y, ec->size, 0)) {
        ak_set_curve_form_e(ec);
        ak_point_set_as_unit(ep, ec);
		return;
	}
	ak_mpznmax a_,b_;
    ak_mpzn_mul_montgomery(a_, wp->z, ec->edw.t, ec->p, ec->n, ec->size);
    ak_mpzn_sub(a_, ec->p, a_, ec->size);
	ak_mpzn_add_montgomery(a_, wp->x, a_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(b_, wp->z, ec->edw.s, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(ep->x, b_, a_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->z, ep->x, wp->y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->x, a_, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_sub(ep->y, ec->p, b_, ec->size);
    ak_mpzn_add_montgomery(ep->y, ep->y, a_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->y, ep->y, wp->y, ec->p, ec->n, ec->size);
    ak_set_curve_form_e(ec);
}
void ak_epoint_triple(ak_point ep, ak_curve ec){
    ak_mpznmax a,b,c,d,e,f,g;
    ak_mpzn_mul_montgomery(a, ep->y, ep->y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, ep->x, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(b, b, ec->a, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(c, a, b, ec->p, ec->size);
    ak_mpzn_mul_montgomery(d, ep->z, ep->z, ec->p, ec->n, ec->size);
    ak_mpzn_lshift_montgomery(d,d,ec->p, ec->size);
    ak_mpzn_sub(e,ec->p,c,ec->size);
    ak_mpzn_add_montgomery(d, d, e, ec->p, ec->size);
    ak_mpzn_lshift_montgomery(d,d,ec->p, ec->size); //b=2(2*Z^2...)
    ak_mpzn_mul_montgomery(f, b, d, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(g, a, d, ec->p, ec->n, ec->size);
    ak_mpzn_sub(e,ec->p,b,ec->size);
    ak_mpzn_add_montgomery(e, a, e, ec->p, ec->size);
    ak_mpzn_mul_montgomery(e, e, c, ec->p, ec->n, ec->size);
    ak_mpzn_sub(a,ec->p,g,ec->size);
    ak_mpzn_add_montgomery(a, e, a, ec->p, ec->size);
    ak_mpzn_add_montgomery(b, e, f, ec->p, ec->size);
    ak_mpzn_add_montgomery(c, g, e, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->x, ep->x, c, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->x, ep->x, a, ec->p, ec->n, ec->size);
    ak_mpzn_sub(c,ec->p,e,ec->size);
    ak_mpzn_add_montgomery(c, c, f, ec->p, ec->size);
    ak_mpzn_mul_montgomery(ep->y, ep->y, c, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->y, ep->y, b, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->z, ep->z, a, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->z, ep->z, b, ec->p, ec->n, ec->size);
}
void ak_wpoint_triple(ak_point wp, ak_curve ec){
    struct point a;
    ak_point wp_=&a;
    ak_point_set_point(wp_,wp,ec);
    ak_wpoint_double((ak_wpoint)wp_,(ak_wcurve)ec);
    ak_wpoint_add((ak_wpoint)wp,(ak_wpoint)wp_,(ak_wcurve)ec);
}
void ak_point_triple(ak_point ep, ak_curve ec){
    if (ec->cform==0) ak_wpoint_triple(ep,ec);
    if (ec->cform==1)   ak_epoint_triple(ep,ec);
}
void ak_epoint_quintuple(ak_point ep, ak_curve ec){
    ak_mpznmax x5,y5,z5,t1,t2,t_;
    ak_mpzn_mul_montgomery(z5, ep->y, ep->y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(y5, ep->x, ep->x, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, ec->a, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(x5, z5, y5, ec->p, ec->size);
    ak_mpzn_mul_montgomery(t1, ep->z, ep->z, ec->p, ec->n, ec->size);
    ak_mpzn_lshift_montgomery(t1,t1,ec->p, ec->size);
    ak_mpzn_sub(t1,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(t1, t1, x5, ec->p, ec->size);
    ak_mpzn_lshift_montgomery(z5,z5,ec->p, ec->size);
    ak_mpzn_lshift_montgomery(y5,y5,ec->p, ec->size);
    ak_mpzn_mul_montgomery(z5, z5, t1, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(t1, t1, y5, ec->p, ec->n, ec->size);
    ak_mpzn_sub(y5,ec->p,y5,ec->size);
    ak_mpzn_add_montgomery(y5, y5, x5, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, x5, y5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(x5, y5, z5, ec->p, ec->size);
    ak_mpzn_sub(t2,ec->p,z5,ec->size);
    ak_mpzn_add_montgomery(t2, y5, t2, ec->p, ec->size);
    ak_mpzn_mul_montgomery(x5, x5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_sub(t2,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(t2, y5, t2, ec->p, ec->size);
    ak_mpzn_add_montgomery(y5, y5, t1, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(t2, t2, t1, ec->p, ec->size);
    ak_mpzn_mul_montgomery(t1, t1, x5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(z5, z5, y5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(x5, x5, t2, ec->p, ec->n, ec->size);
    ak_mpzn_sub(x5,ec->p,x5,ec->size);
    ak_mpzn_mul_montgomery(y5, t2, y5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(t2, x5, z5, ec->p, ec->size);
    ak_mpzn_sub(t_,ec->p,z5,ec->size);
    ak_mpzn_add_montgomery(x5, x5, t_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(x5, t2, x5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->x, ep->x, x5, ec->p, ec->n, ec->size);
    ak_mpzn_add_montgomery(z5, y5, t1, ec->p, ec->size);
    ak_mpzn_sub(t_,ec->p,t1,ec->size);
    ak_mpzn_add_montgomery(y5, y5, t_, ec->p, ec->size);
    ak_mpzn_mul_montgomery(y5, y5, z5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->y, y5, ep->y, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(z5, t2, z5, ec->p, ec->n, ec->size);
    ak_mpzn_mul_montgomery(ep->z, ep->z, z5, ec->p, ec->n, ec->size);
}
void ak_wpoint_quintuple(ak_point wp, ak_curve ec){
    struct point a;
    ak_point wp_=&a;
    ak_point_set_point(wp_,wp,ec);
    ak_wpoint_double((ak_wpoint)wp_,(ak_wcurve)ec);
    ak_wpoint_double((ak_wpoint)wp_,(ak_wcurve)ec);
    ak_wpoint_add((ak_wpoint)wp,(ak_wpoint)wp_,(ak_wcurve)ec);}
void ak_point_quintuple(ak_point ep, ak_curve ec){
    if (ec->cform==0) ak_wpoint_quintuple(ep,ec);
    if (ec->cform==1)   ak_epoint_quintuple(ep,ec);
}

void ak_point_pow_binary( ak_point wq, ak_point wp, ak_uint64 *k, size_t size,
                          ak_curve ec ){
    ak_uint64 uk;
    long long int i, j;
    struct point Q, R;
    ak_point_set_as_unit(&Q, ec);
    ak_point_set_point( &R, wp, ec );
    for( i = size-1; i >= 0; i-- ) {
        uk = k[i];
        for (j = 0; j < 64; j++) {
            ak_point_double( &Q, ec );
            if (uk & 0x8000000000000000LL) ak_point_add( &Q, &R, ec );
            uk <<= 1;
        }
    }
    ak_point_set_point( wq, &Q, ec );
}
ak_uint32 ak_mpzn_rem_uint32( ak_uint64 *x, const size_t size, ak_uint32 p )
{
    size_t i;
    ak_uint64 t, r1, r = 1, sum = x[0]%p;
    if( !p ) return ak_error_message( ak_error_invalid_value, __func__, "divide by zero" );
    r1 = 9223372036854775808ull % p; r1 = ( 2*r1 )%p;
    for( i = 1; i < size; i++ ) {
        r *= r1; r %= p;
        t = x[i]%p; t *= r; t %= p; sum += t;
    }
    return sum%p;
}

ak_int32 ak_mods(ak_uint64* n,ak_uint32 l,ak_uint32 w,size_t size){
    ak_uint32 a=(pow(l,w-1));
    ak_uint32 p=a*l;
    ak_uint32 res=ak_mpzn_rem_uint32(n, size, p); //result = n mod l^w
    if(res>=a) return res-p;
    return res;
}
void ak_mpzn_div2(ak_uint64* n, size_t size){
    for(unsigned i=0; i<size-1; i++){
        n[i]>>=1;
        if(n[i+1]&1) n[i]|=0x8000000000000000LL;
    }
    n[size-1]>>=1;
}
ak_int32 ak_n_to_NAF2(ak_uint64 *k, ak_int32* res, size_t size){
    /*
     * Atomicity Improvement for Elliptic Curve Scalar Multiplication
     */
    ak_int32 i=0;
    ak_mpznmax E,one=ak_mpznmax_one;
    ak_mpzn_set(E,k,size);
    while (!ak_mpzn_cmp_ui(E,size,0)){
        if(E[0]&1) {
            res[i]=2-(E[0]&3);
            if (res[i]==1) E[0]-=1;
            else (ak_mpzn_add(E,E,one, size));
        }
        else res[i]=0;
        ak_mpzn_div2(E,size);
        i++;
    }
    return i;
}
void ak_invert_p(ak_point resp,ak_point ep, ak_curve ec){
    if(ec->cform==0){
        ak_mpzn_sub(resp->y,ec->p,ep->y,ec->size);
        memcpy(resp->x, ep->x, ec->size * sizeof(ak_uint64));
        memcpy(resp->z, ep->z, ec->size * sizeof(ak_uint64));
    }
    if(ec->cform==1){
        ak_mpzn_sub(resp->x,ec->p,ep->x,ec->size);
        memcpy(resp->y, ep->y, ec->size * sizeof(ak_uint64));
        memcpy(resp->z, ep->z, ec->size * sizeof(ak_uint64));
    }
}
void ak_point_pow_NAF(ak_point wq, ak_point wp, ak_int32 *kNAF,ak_int32 i,
                      ak_curve ec){
    struct point Q, R,R_;
    ak_point_set_as_unit( &Q, ec );
    ak_point_set_point( &R, wp, ec );
    ak_point_set_point( &R_, wp, ec );
    ak_invert_p(&R_,&R,ec);
    for (i=i-1; i>=0; i--) {
        ak_point_double(&Q,ec);
        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_point_add(&Q,&R,ec);
        }
        else if(nj<0){
            ak_point_add(&Q,&R_,ec);
        }
    }
    ak_point_set_point( wq, &Q, ec );
}
void ak_point_pow_NAF_wv(ak_point wq, ak_point wp, ak_uint64 *k, size_t size,
                       ak_curve ec)
{
    ak_int32 kNAF[512];
    ak_int32 i=ak_n_to_NAF2(k,kNAF,size);
    ak_point_pow_NAF(wq, wp, kNAF,i,ec);
}

ak_int32 ak_n_to_NAF_powof2(ak_uint64 *k,ak_int32 w, ak_int32* res, size_t size){
    ak_int32 i=0,mods=pow(2,w);
    ak_mpznmax E,one=ak_mpznmax_one;
    ak_mpzn_set(E,k,size);
    while (!ak_mpzn_cmp_ui(E,size,0)){
        if(E[0]&1){
            res[i]=(E[0]&(mods-1));
            if(res[i]>=(mods)>>1) res[i]-=(mods);
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }
        else res[i]=0;
        ak_mpzn_div2(E,size);
        i++;
    }
    return i;
}
void ak_point_pow_NAF_powof2( ak_point wq, ak_point wp, ak_int32 *kNAF, ak_int32 i,ak_uint32 w,
                              ak_curve ec){
    ak_mpznmax pows;
    struct point Q, R_[32],inv;
    ak_point_set_as_unit(&Q, ec );
    ak_point_set_point(R_, wp, ec );
    ak_point_set_point(R_+1, wp, ec );
    ak_point_set_point(R_+2, wp, ec );
    ak_point_triple(R_+1, ec);
    if(w>3){
        ak_point_quintuple(R_+2, ec);
        for (ak_uint32 j=3;j<w;j++){
            ak_mpzn_set_ui(pows,ec->size,j*2+1);
            ak_point_pow_binary(R_+j,R_,pows,ec->size,ec);
        }
    }
    for (i=i-1; i>=0; i--) {
        ak_point_double(&Q,ec);
        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_point_add(&Q,R_+((nj-1)>>1),ec);
        }
        else if(nj<0){
            ak_invert_p(&inv,R_+((-nj-1)>>1),ec);
            ak_point_add(&Q,&inv,ec);
        }
    }
    ak_point_set_point( wq, &Q, ec );
}

void ak_point_pow_NAF_powof2_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size,ak_uint32 w,
                       ak_curve ec)
{/*
Fixed-Base Comb with Window-Non-Adjacent Form (NAF)
Method for Scalar Multiplication
*/
    ak_int32 kNAF[512];
    ak_int32 i=ak_n_to_NAF_powof2(k,w,kNAF,size);
    ak_point_pow_NAF_powof2(wq, wp, kNAF,i,w, ec);

}


ak_int32 ak_n_to_NAF_L_w(ak_uint64 *k, ak_int32* res, ak_int32 L, ak_int32 w, size_t size,ak_uint64 *p,ak_uint64 *r2,ak_uint64 n){
    ak_int32 i=0,L_w1=pow(L,w-1);
    ak_int32 L_w=L_w1*L;
    ak_mpznmax E,one,L_inv,u=ak_mpznmax_one,z;
    ak_mpzn_set_ui(z, size, 2);
    ak_mpzn_sub(z, p, z, size);
    ak_mpzn_set_ui(L_inv,size,L);
    ak_mpzn_mul_montgomery(L_inv,L_inv,r2,p,n,size);
    ak_mpzn_modpow_montgomery(L_inv, L_inv, z, p, n, size);

    ak_mpzn_set(E,k,size);


    while (ak_mpzn_cmp(E,u,size)>=0){

        if(ak_mpzn_rem_uint32(E,size,L)!=0){
            res[i]=(ak_mpzn_rem_uint32(E,size,L_w));
            if(res[i]>=L_w1) res[i]-=(L_w);
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                //if(!ak_mpzn_cmp(E,one,size)) break;
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }
        else res[i]=0;
        ak_mpzn_mul_montgomery(E,E,r2,p,n,size);
        ak_mpzn_mul_montgomery(E,E,L_inv,p,n,size);
        ak_mpzn_mul_montgomery(E,E,u,p,n,size);
        i++;
    }

    return i;
}
void ak_point_pow_NAF_powofL( ak_point wq, ak_point wp, ak_int32 *kNAF, ak_int32 i, ak_uint32 l,ak_uint32 w,
                              ak_curve ec){
    ak_mpznmax pows;
    ak_int32 LL_=pow(l,w-1);
    ak_int32 len=LL_*(l-1);
    struct point Q, R_[12500],inv;
    ak_point_set_as_unit(&Q, ec );
    ak_point_set_point(R_, wp, ec );
    ak_int32 cr=1;
    for (ak_int32 j=2;j<len;j++){
        if(j%l){ ak_mpzn_set_ui(pows,ec->size,j);
            ak_point_pow_binary(R_+cr,R_,pows,ec->size,ec);
            cr++;
        }
    }

    for (i=i-1; i>=0; i--) {
        if(l==2) ak_point_double(&Q,ec);
        else{
            if(l==3) ak_point_triple(&Q,ec);
            else if(l==5) ak_point_quintuple(&Q,ec);
            else{
                ak_mpzn_set_ui(pows,ec->size,l);
                ak_point_pow_binary(&Q,&Q,pows,ec->size,ec);
            }
        }

        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_point_add(&Q,R_+(nj-1-(nj/l)),ec);
        }
        else if(nj<0){
            ak_invert_p(&inv,R_+((-nj)-1-(-nj)/l),ec);
            ak_point_add(&Q,&inv,ec);
        }
    }

    ak_point_set_point( wq, &Q, ec );
}
void ak_point_pow_NAF_powofL_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size, ak_uint32 l,ak_uint32 w,
                              ak_curve ec){
    ak_int32 kNAF[512];
    ak_int32 i=ak_n_to_NAF_L_w(k,kNAF,l,w,size,ec->p,ec->r2,ec->n);
    ak_point_pow_NAF_powofL(wq,wp,kNAF,i,l,w,ec);
}

ak_int32 ak_n_to_Ext_wmb_NAF(ak_uint64 *k, ak_int32* res, ak_int32* bases,
                        ak_uint32* a, ak_uint32* w, ak_int8 len,size_t size,ak_uint64 *p,ak_uint64 *r2,ak_uint64 n){

    ak_int32 i=0;
    ak_int32 A=1;
    ak_mpznmax E,one,L_inv[10],u=ak_mpznmax_one,z;
    for(int j=0;j<len;j++) A*=pow(a[j],w[j]);

    ak_mpzn_set_ui(z, size, 2);
    ak_mpzn_sub(z, p, z, size);
    for(int j=0;j<len;j++){
        ak_mpzn_set_ui(L_inv[j],size,a[j]);
        ak_mpzn_mul_montgomery(L_inv[j],L_inv[j],r2,p,n,size);
        ak_mpzn_modpow_montgomery(L_inv[j], L_inv[j], z, p, n,size);
    }

    ak_mpzn_set(E,k,size);


    while (ak_mpzn_cmp(E,u,size)>=0){
        int j=0;
        for(;j<len;j++){
            if(ak_mpzn_rem_uint32(E,size,a[j])==0) {res[i]=0; break;}
        }
        if(j==len) {
            res[i]=ak_mpzn_rem_uint32(E,size,A);
            if(res[i]>=A/2) res[i]-=A;
            if (res[i]>0){
                ak_mpzn_set_ui(one,size,res[i]);
                ak_mpzn_sub(E,E,one, size);
            }
            else{
                ak_mpzn_set_ui(one,size,-res[i]);
                ak_mpzn_add(E,E,one, size);
            }
        }

        for(j=0;j<len;j++){
            if(ak_mpzn_rem_uint32(E,size,a[j])==0){
                ak_mpzn_mul_montgomery(E,E,r2,p,n,size);
                ak_mpzn_mul_montgomery(E,E,L_inv[j],p,n,size);
                ak_mpzn_mul_montgomery(E,E,u,p,n,size);
                bases[i]=a[j];
                break;
            }
        }
        i++;
    }

    return i;
}
void ak_point_pow_NAF_mbw( ak_point wq, ak_point wp, ak_int32 *kNAF, ak_int32 i,ak_int32 *bases, ak_uint32* l,ak_uint32* w, ak_uint8 len,
                             ak_curve ec){
    ak_mpznmax pows;
    ak_int32 A=1;
    for(int j=0;j<len;j++) A*=pow(l[j],w[j]);
    struct point Q, R_[25500],inv;
    ak_point_set_as_unit(&Q, ec );
    ak_point_set_point(R_, wp, ec );
    for (ak_int32 j=2;j<A/2;j++){
        ak_mpzn_set_ui(pows,ec->size,j);
        ak_point_pow_binary(R_+j-1,R_,pows,ec->size,ec);
    }

    for (i=i-1; i>=0; i--) {

        switch (bases[i]) {
        case 2:
            ak_point_double(&Q,ec); break;
        case 3:
            ak_point_triple(&Q,ec); break;
        case 5:
            ak_point_quintuple(&Q,ec); break;
        default:
            ak_mpzn_set_ui(pows,ec->size,bases[i]);
            ak_point_pow_binary(&Q,&Q,pows,ec->size,ec);
            break;
        }

        ak_int32 nj=kNAF[i];
        if (nj>0) {
            ak_point_add(&Q,R_+(nj-1),ec);
        }
        else if(nj<0){
            ak_invert_p(&inv,R_+(-nj-1),ec);
            ak_point_add(&Q,&inv,ec);
        }
    }

    ak_point_set_point( wq, &Q, ec );

}
void ak_point_pow_NAF_mbw_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size, ak_uint32* l,ak_uint32* w, ak_uint8 len,
                                   ak_curve ec){
    ak_int32 kNAF[512],bases[512];
    ak_int32 i=ak_n_to_Ext_wmb_NAF(k,kNAF,bases,l,w,len,size,ec->p,ec->r2,ec->n);
    ak_point_pow_NAF_mbw( wq,wp,kNAF,i,bases,l,w,len,ec);
}
