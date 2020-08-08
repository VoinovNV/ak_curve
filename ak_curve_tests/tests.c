#include "curve.h"
#include <assert.h>
#include "tests.h"
#include "ak_curves.h"
#include <ak_parameters.h>
#include <math.h>
#include <unistd.h>
int cmp_point(ak_point ep1, ak_point ep2, ak_curve ec) {
    if (ak_mpzn_cmp(ep1->x, ep2->x, ec->size)) return 1;
    if (ak_mpzn_cmp(ep1->y, ep2->y, ec->size)) return 2;
    if (ak_mpzn_cmp(ep1->z, ep2->z, ec->size)) return 3;
    return 0;
}
void test_add_(ak_curve ec) {
    struct point ep1, ep2, ep3, ep4;
    ak_point point1 = &ep1, point2 = &ep2, point3 = &ep3, point4 = &ep4;

    ak_point_set_point(point1, &ec->point, ec);
    ak_point_set_as_unit(point2, ec);
    ak_point_set_as_unit(point3, ec);
    ak_point_set_point(point4, point1, ec);

    ak_point_add(point1, point2, ec);
    ak_point_reduce(point1, ec);

    assert(ak_point_is_ok(point1, ec) && "test add : P + 0 : Point is not on a curve!\n");
    assert(!cmp_point(point1, point4, ec) && "test add : P + 0 : Something is wrong!\n");

    ak_point_add(point2, point2, ec);
    ak_point_reduce(point2, ec);

    assert(ak_point_is_ok(point2, ec) && "test add : 0 + 0 : Point is not on a curve!\n");
    assert(!cmp_point(point2, point3, ec) && "test add : 0 + 0 : Something is wrong!\n");

}
void test_add(){
    struct curve gost_e = gost_3410_2012_512_paramSetC,gost_w=gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost_e,wc=&gost_w;
    test_add_(ec);
    test_add_(wc);
}
void test_double_(ak_curve ec) {
    struct point ep1, ep2;
    ak_point point1 = &ep1, point2 = &ep2;

    ak_point_set_as_unit(point1, ec);
    ak_point_set_as_unit(point2, ec);

    ak_point_double(point1, ec);
    ak_point_reduce(point1, ec);

    assert(ak_point_is_ok(point1, ec) && "test double : 2*0 = 0 : Point is not on a curve!\n");
    assert(!cmp_point(point1, point2, ec) && "test double : 2*0 = 0 : Something is wrong!\n");

}
void test_double(){
    struct curve gost_e = gost_3410_2012_512_paramSetC,gost_w=gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost_e,wc=&gost_w;
    test_double_(ec);
    test_double_(wc);
}
void test_pow(void (*f)( ak_point , ak_point, ak_uint64 *, size_t,
                                       ak_curve, ak_uint64, ak_uint64 ),ak_uint64 l, ak_uint64 w) {
    struct curve gost = gost_3410_2012_512_paramSetC, gost_w = gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost;

    for(int qwe=0;qwe<2;qwe++){
    struct point ep1, ep2, ep3, ep4;
    ak_point point1 = &ep1, point2 = &ep2, point3 = &ep3, point4 = &ep4;

    ak_point_set_point(point1, &ec->point, ec);
    ak_point_set_as_unit(point2, ec);
    ak_point_set_point(point3, point1, ec);
    ak_point_set_point(point4, point1, ec);

    //ak_point_pow(point1, point1, ec->q, ec->size, ec);
    f(point1, point1, ec->q, ec->size, ec,l,w);
    ak_point_reduce(point1, ec);

    assert(ak_point_is_ok(point1, ec) && "test pow : qP=0 : Point is not on a curve!\n");
    assert(!cmp_point(point1, point2, ec) && "test pow : qP=0 : Something is wrong!");

	ak_mpzn512 one = ak_mpzn512_one, two, sum;
	ak_mpzn_add_montgomery(one, one, ec->q, ec->p, ec->size);
    //ak_point_pow(point3, point4, one, ec->size, ec);
    f(point3, point4, one, ec->size, ec,l,w);
    ak_point_reduce(point3, ec);
    assert(ak_point_is_ok(point3, ec) && "test pow : (q+1)P=P : Point is not on a curve!\n");
    assert(!cmp_point(point3, point4, ec) && "test pow : (q+1)P=P : Something is wrong!");

	struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 3; i++) {
        ak_point_set_point(point3, &ec->point, ec);
        ak_point_set_point(point4, &ec->point, ec);
		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        //ak_point_pow(point3, point3, one, ec->size, ec);
        f(point3, point3, one, ec->size, ec,l,w);
        ak_point_double(point3, ec);
        ak_point_reduce(point3, ec);
        ak_point_double(point4, ec);
        //ak_point_pow(point4, point4, one, ec->size, ec);
        f(point4, point4, one, ec->size, ec,l,w);
        ak_point_reduce(point4, ec);
        assert(ak_point_is_ok(point3, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 1 is not on a curve!\n");
        assert(ak_point_is_ok(point4, ec) && "test pow : 2((r1)P)=r1((2)P) : Point 2 is not on a curve!\n");
        assert(!cmp_point(point3, point4, ec) && "test pow : 2((r1)P)=r1((2)P) : Something is wrong!");

		ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
		ak_mpzn_set_random_modulo(two, ec->q, ec->size, &r);
        ak_point_set_point(point3, &ec->point, ec);
        ak_point_set_point(point4, &ec->point, ec);
		ak_mpzn_add_montgomery(sum, one, two, ec->p, ec->size);


        //ak_point_pow(point4, point4, sum, ec->size, ec);
        f(point4, point4, sum, ec->size, ec,l,w);
        ak_point_reduce(point4, ec);
        //ak_point_pow(point1, point3, one, ec->size, ec);
        //ak_point_pow(point2, point3, two, ec->size, ec);
        f(point1, point3, one, ec->size, ec,l,w);
        f(point2, point3, two, ec->size, ec,l,w);
        ak_point_add(point1, point2, ec);
        ak_point_reduce(point1, ec);

        assert(ak_point_is_ok(point1, ec) && "test pow : (r1+r2)P=r1P+r2P : Point 1 is not on a curve!\n");
        assert(ak_point_is_ok(point4, ec) && "test pow : (r1+r2)P=r1P+r2P : Point 2 is not on a curve!\n");
        assert(!cmp_point(point1, point4, ec) && "test pow : (r1+r2)P=r1P+r2P : Something is wrong!");
	}
	ak_random_context_destroy(&r);
    ec = &gost_w;
    }
}
void test_form_changing(){
    struct curve e_gost = gost_3410_2012_512_paramSetC; ak_curve ec = &e_gost;
    struct curve w_gost = gost_3410_2012_512_paramSetC_w; ak_curve wc = &w_gost;
    struct point ep1, ep2;
    ak_point point1 = &ep1, point2 = &ep2;
    struct point wp1, wp2;
    ak_point wpoint1 = &wp1, wpoint2 = &wp2;

    ak_point_set_point(point1, &ec->point, ec);
    ak_point_set_as_unit(point2, ec);

    ak_mpzn512 one=ak_mpzn512_one;

    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 50; i++) {
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_point_set_point(point1, &ec->point, ec);
        ak_point_set_point(wpoint1, &wc->point, wc);

        ak_point_pow(point1, point1, one, ec->size, ec);
        ak_point_pow(wpoint1, wpoint1, one, wc->size, wc);

        ak_epoint_to_wpoint(wpoint2,point1, wc);
        ak_wpoint_to_epoint(point2,wpoint1, ec);
        ak_point_reduce(point2, wc);
        ak_point_reduce(wpoint2, ec);
        ak_point_reduce(point1, wc);
        ak_point_reduce(wpoint1, ec);
        assert(!ak_mpzn_cmp(wpoint1->x, wpoint2->x, ec->size) && "test form : wpoint->x : Something is wrong!");
        assert(!ak_mpzn_cmp(wpoint1->y, wpoint2->y, ec->size) && "test form : wpoint->y : Something is wrong!");
        assert(!ak_mpzn_cmp(point1->x, point2->x, ec->size) && "test form : point->x : Something is wrong!");
        assert(!ak_mpzn_cmp(point1->y, point2->y, ec->size) && "test form : point->y : Something is wrong!");
    }
    ak_random_context_destroy(&r);
}

void test_triple_(ak_curve ec){

    struct point ep1, ep2;
    ak_point point1 = &ep1, point2 = &ep2;
    ak_mpzn512 one = ak_mpzn512_one;

    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 50; i++) {
        ak_point_set(point1, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_point_pow(point1, point1, one, ec->size, ec);
        ak_point_set_point(point2,point1,ec);
        ak_point_triple(point1, ec);
        ak_mpzn_set_ui(one,ak_mpzn512_size,3);
        ak_point_pow(point2, point2, one, ec->size, ec);
        ak_point_reduce(point1, ec);
        ak_point_reduce(point2, ec);
        assert(ak_point_is_ok(point1, ec) && "test trip : Point 1 is not on a curve!\n");
        assert(ak_point_is_ok(point2, ec) && "test trip : Point 2 is not on a curve!\n");
        assert(!cmp_point(point1, point2, ec) && "test trip : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}
void test_triple(){
    struct curve gost_e = gost_3410_2012_512_paramSetC,gost_w=gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost_e,wc=&gost_w;
    test_triple_(ec);
    test_triple_(wc);
}
void test_quintuple_(ak_curve ec){

    struct point ep1, ep2;
    ak_point point1 = &ep1, point2 = &ep2;
    ak_mpzn512 one = ak_mpzn512_one;

    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 1; i++) {
        ak_point_set(point1, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_point_pow(point1, point1, one, ec->size, ec);
        ak_point_set_point(point2,point1,ec);
        ak_point_quintuple(point1, ec);
        ak_mpzn_set_ui(one,ak_mpzn512_size,5);
        ak_point_pow(point2, point2, one, ec->size, ec);
        ak_point_reduce(point1, ec);
        ak_point_reduce(point2, ec);

        assert(ak_point_is_ok(point1, ec) && "test trip : Point 1 is not on a curve!\n");
        assert(ak_point_is_ok(point2, ec) && "test trip : Point 2 is not on a curve!\n");
        assert(!cmp_point(point1, point2, ec) && "test trip : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}
void test_quintuple(){
    struct curve gost_e = gost_3410_2012_512_paramSetC,gost_w=gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost_e,wc=&gost_w;
    test_quintuple_(ec);
    test_quintuple_(wc);
}
void M_pow_mont( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow(wq,ep,k,size,ec);
}
void test_pow_montgomery(){
    test_pow(M_pow_mont,0,0);
}
void M_pow_bin( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_binary(wq,ep,k,size,ec);
}
void test_pow_bin(){
    test_pow(M_pow_bin,0,0);
}
void M_pow_NAF( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_wv(wq,ep,k,size,ec);
}

void test_pow_NAF(){
    test_pow(M_pow_NAF,0,0);
}
void M_pow_NAF_2_w( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_powof2_wv(wq,ep,k,size,w,ec);
}
void test_pow_NAF_2_w(){
    for(int i=2;i<5;i++)test_pow(M_pow_NAF_2_w,0,i);
}


void M_pow_NAF_l_w( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_powofL_wv(wq,ep,k,size,l,w,ec);
}
void test_pow_NAF_l_w(){
    for(int i=2;i<5;i++) test_pow(M_pow_NAF_l_w,3,i);
    for(int i=2;i<5;i++) test_pow(M_pow_NAF_l_w,5,i);
}


void test_pow_NAF_mbw(){
    struct curve gost = gost_3410_2012_512_paramSetC;
    ak_curve ec = &gost;
    struct point ep1, ep2;
    ak_point point1 = &ep1, point2 = &ep2;
    const int ss=2;
    ak_uint32 l[ss]={2,3};
    ak_uint32 w[ss]={1,2};
    ak_point_set_point(point1, &ec->point, ec);

    ak_mpzn512 one;
    struct random r;
    ak_random_context_create_lcg(&r);
    for (int i = 0; i <= 25; i++) {
        ak_point_set_point(point1, &ec->point, ec);
        ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
        ak_point_pow(point2, point1, one, ec->size, ec);
        ak_point_pow_NAF_mbw_wv(point1, point1, one, ec->size,l,w,ss,ec);
        ak_point_reduce(point1, ec);
        ak_point_reduce(point2, ec);
        assert(!cmp_point(point1, point2, ec) && "MBW pow : Something is wrong!");

    }
    ak_random_context_destroy(&r);
}
