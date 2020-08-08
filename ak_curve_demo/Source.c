#include "curve.h"
#include <math.h>
#include <x86intrin.h>
#define NUM_OF_ITER 99
unsigned long long GET_TIME_ALL(void f(ak_point, ak_point, ak_uint64*, size_t, ak_curve, ak_uint64, ak_uint64),ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w, int cnt) { \
    unsigned long long a=_rdtsc();
    for(int i=0;i<cnt;i++) f(wq,ep,k,size,ec,l,w);
    unsigned long long b=_rdtsc();
    return b-a;
}
void M_pow_mont( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow(wq,ep,k,size,ec);
}
void M_pow_bin( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_binary(wq,ep,k,size,ec);
}
void M_pow_NAF( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_wv(wq,ep,k,size,ec);
}
void M_pow_NAF_2_w( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_powof2_wv(wq,ep,k,size,w,ec);
}
void M_pow_NAF_l_w( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    ak_point_pow_NAF_powofL_wv(wq,ep,k,size,l,w,ec);
}
void M_pow_mont_ch_f( ak_point wq, ak_point ep, ak_uint64 *k, size_t size, ak_curve ec, ak_uint64 l, ak_uint64 w){
    if(ec->cform==0) ak_wpoint_to_epoint(wq,ep,ec);
    else ak_epoint_to_wpoint(wq,ep,ec);
    ak_point_pow(wq,wq,k,size,ec);
    if(ec->cform==0) ak_wpoint_to_epoint(ep,wq,ec);
    else ak_epoint_to_wpoint(ep,wq,ec);

}


void all_time(ak_curve ec,ak_mpzn512 one){
    const int num_func=6;
    void (*f[num_func])(ak_point, ak_point, ak_uint64*, size_t, ak_curve, ak_uint64, ak_uint64)={
        M_pow_mont,
        M_pow_bin,
        M_pow_NAF,
        M_pow_NAF_2_w,
        M_pow_NAF_l_w,
        M_pow_mont_ch_f
            };
    char* str[num_func]={
        "Лесенка Монтгомери",
        "Бинарный алгоритм",
        "С основанием 2",
        "С основанием 2 и окном",
        "С основание L и окном",
        "С переходом к другой форме",
        };
    ak_uint64 l[num_func][4]={
        {1,0},{1,0},{1,0},{1,0},{3,5,3,2},{1,0}
    };
    ak_uint64 w[num_func][4]={
        {1,0},{1,0},{1,0},{3,2,3,5},{3,2,3,5},{1,0}
    };

    struct point un_p;
    ak_point ep=&un_p;
    for(int i=0;i<num_func;i++)
        for(unsigned long long j=1;j<=l[i][0];j++)
            for(unsigned long long k=1;k<=w[i][0];k++)
            {
                ak_point_set(ep,ec);
                for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1))
                    printf("%s: (основание: %llu окно:%llu);%d ; %llu\n",str[i],l[i][j],w[i][k],cnt,GET_TIME_ALL(f[i],ep,ep,one,ak_mpzn512_size,ec,l[i][j],w[i][k],cnt));
            }
}

void all_time_mbw (ak_curve ec,ak_mpzn512 one){
    struct point un_p;
    ak_point ep=&un_p;
    for(int i=1;i<4;i++){
        ak_uint32 l[2]={2,3};
        ak_uint32 w[2]={i,i};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_mbw_wv(ep,ep,one,ak_mpzn512_size,l,w,2,ec);
            unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

    }
    for(int i=1;i<4;i++){
        ak_uint32 l[3]={2,3,5};
        ak_uint32 w[3]={i,i,i};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_mbw_wv(ep,ep,one,ak_mpzn512_size,l,w,3,ec);
            unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3,5 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

    }
}
void time_with_point_eval(ak_curve ec,ak_mpzn512 one){
    struct point un_p;
    ak_point ep=&un_p;

    {
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            ak_int32 kNAF[512];
            ak_int32 i=ak_n_to_NAF2(one,kNAF,ak_mpzn512_size);
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF(ep, ep, kNAF,i,ec);
            unsigned long long b=_rdtsc();
            printf("С основанием 2 и окном: (основание: %d);%d ; %llu\n",2,cnt,b-a);
        }
    }
    {

        ak_uint64 w[3]={2,3,5};
        for(unsigned long long k=0;k<3;k++)
        {
                ak_point_set(ep,ec);
                for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
                    ak_int32 kNAF[512];
                    ak_int32 i=ak_n_to_NAF_powof2(one,w[k],kNAF,ak_mpzn512_size);
                    unsigned long long a=_rdtsc();
                    for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_powof2(ep, ep, kNAF,i,w[k], ec);
                    unsigned long long b=_rdtsc();
                    printf("С основание 2 и окном: (основание: %d окно:%llu);%d ; %llu\n",2,w[k],cnt,b-a);
                }
        }


    }
    {
    ak_uint64 l[3]={2,3,5};
    ak_uint64 w[3]={2,3,5};
    for(unsigned long long j=0;j<3;j++)
        for(unsigned long long k=0;k<3;k++)
            {
                ak_point_set(ep,ec);
                for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
                    ak_int32 kNAF[512];
                    ak_int32 i=ak_n_to_NAF_L_w(one,kNAF,l[j],w[k],ak_mpzn512_size,ec->p,ec->r2,ec->n);

                unsigned long long a=_rdtsc();
                for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_powofL(ep,ep,kNAF,i,l[j],w[k],ec);
                unsigned long long b=_rdtsc();
                printf("С основание L и окном: (основание: %llu окно:%llu);%d ; %llu\n",l[j],w[k],cnt,b-a);
            }
        }
    }
    for(int d=1;d<4;d++){
        ak_uint32 l[2]={2,3};
        ak_uint32 w[2]={d,d};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            ak_int32 kNAF[512],bases[512];
            ak_int32 i=ak_n_to_Ext_wmb_NAF(one,kNAF,bases,l,w,2,ak_mpzn512_size,ec->p,ec->r2,ec->n);
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_mbw( ep,ep,kNAF,i,bases,l,w,2,ec);
            unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

    }
    for(int d=1;d<4;d++){
        ak_uint32 l[3]={2,3,5};
        ak_uint32 w[3]={d,d,d};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            ak_int32 kNAF[512],bases[512];
            ak_int32 i=ak_n_to_Ext_wmb_NAF(one,kNAF,bases,l,w,3,ak_mpzn512_size,ec->p,ec->r2,ec->n);
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_mbw( ep,ep,kNAF,i,bases,l,w,3,ec);
            unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3,5 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

    }
}
void time_only_pow(ak_curve ec,ak_mpzn512 one){
    struct point un_p;
    ak_point ep=&un_p;

    {
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            ak_int32 kNAF[512];
            ak_int32 i=ak_n_to_NAF2(one,kNAF,ak_mpzn512_size);
            struct point Q, R,R_;
            ak_point_set_as_unit( &Q, ec );
            ak_point_set_point( &R, ep, ec );
            ak_point_set_point( &R_, ep, ec );
            ak_invert_p(&R_,&R,ec);

            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) {for (i=i-1; i>=0; i--) {
                    ak_point_double(&Q,ec);
                    ak_int32 nj=kNAF[i];
                    if (nj>0) {
                        ak_point_add(&Q,&R,ec);
                    }
                    else if(nj<0){
                        ak_point_add(&Q,&R_,ec);
                    }
                }
                ak_point_set_point( ep, &Q, ec );
            }
            unsigned long long b=_rdtsc();
            printf("С основанием 2;%d ; %llu\n",cnt,b-a);
        }
    }
    {

        ak_uint64 w[3]={2,3,5};
        for(unsigned long long k=0;k<3;k++)
        {
            ak_point_set(ep,ec);
            for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
                ak_int32 kNAF[512];
                ak_int32 i=ak_n_to_NAF_powof2(one,w[k],kNAF,ak_mpzn512_size);
                unsigned long long a=_rdtsc();
                for(int qwe=0;qwe<cnt;qwe++) ak_point_pow_NAF_powof2(ep, ep, kNAF,i,w[k], ec);
                unsigned long long b=_rdtsc();
                printf("С основание 2 и окном: (основание: %d окно:%llu);%d ; %llu\n",2,w[k],cnt,b-a);
            }
        }
    }
    {
        ak_uint64 l[3]={2,3,5};
        ak_uint64 w[3]={2,3,5};
        for(unsigned long long j=0;j<3;j++)
            for(unsigned long long k=0;k<3;k++)
            {
                ak_point_set(ep,ec);
                for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
                    ak_int32 kNAF[10000];
                    ak_int32 i=ak_n_to_NAF_L_w(one,kNAF,l[j],w[k],ak_mpzn512_size,ec->p,ec->r2,ec->n);
                    //
                    ak_mpznmax pows;
                    ak_int32 LL_=pow(l[j],w[k]-1);
                    ak_int32 len=LL_*(l[j]-1);
                    struct point Q, R_[2500],inv;
                    ak_point_set_as_unit(&Q, ec );
                    ak_point_set_point(R_, ep, ec );
                    ak_int32 cr=1;
                    for (ak_int32 d=2;d<len;d++){
                        if(d%l[j]){ ak_mpzn_set_ui(pows,ec->size,d);
                            ak_point_pow_binary(R_+cr,R_,pows,ec->size,ec);
                            cr++;
                        }
                    }

                    //
                    unsigned long long a=_rdtsc();
                    for(int qwe=0;qwe<cnt;qwe++) {for (i=i-1; i>=0; i--) {
                            if(l[j]==2) ak_point_double(&Q,ec);
                            else{
                                if(l[j]==3) ak_point_triple(&Q,ec);
                                else if(l[j]==5) ak_point_quintuple(&Q,ec);
                                else{
                                    ak_mpzn_set_ui(pows,ec->size,l[j]);
                                    ak_point_pow_binary(&Q,&Q,pows,ec->size,ec);
                                }
                            }

                            ak_int32 nj=kNAF[i];
                            if (nj>0) {
                                ak_point_add(&Q,R_+(nj-1-(nj/l[j])),ec);
                            }
                            else if(nj<0){
                                ak_invert_p(&inv,R_+((-nj)-1-(-nj)/l[j]),ec);
                                ak_point_add(&Q,&inv,ec);
                            }
                        }
                        ak_point_set_point( ep, &Q, ec );}
                    unsigned long long b=_rdtsc();
                    printf("С основание L и окном: (основание: %llu окно:%llu);%d ; %llu\n",l[j],w[k],cnt,b-a);
                }
            }
    }
    for(int d=1;d<4;d++){
        ak_uint32 l[2]={2,3};
        ak_uint32 w[2]={d,d};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
        {
                ak_int32 kNAF[512],bases[512];
                ak_int32 i=ak_n_to_Ext_wmb_NAF(one,kNAF,bases,l,w,2,ak_mpzn512_size,ec->p,ec->r2,ec->n);
                //
                ak_mpznmax pows;
                ak_int32 A=1;
                for(int j=0;j<2;j++) A*=pow(l[j],w[j]);
                struct point Q, R_[13500],inv;
                ak_point_set_as_unit(&Q, ec );
                ak_point_set_point(R_, ep, ec );
                for (ak_int32 j=2;j<A/2;j++){
                    ak_mpzn_set_ui(pows,ec->size,j);
                    ak_point_pow_binary(R_+j-1,R_,pows,ec->size,ec);
                }

                //
                unsigned long long a=_rdtsc();
                for(int qwe=0;qwe<cnt;qwe++) {            for (i=i-1; i>=0; i--) {

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

                    ak_point_set_point( ep, &Q, ec );
                }
                unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

        }
    }
    for(int d=1;d<4;d++){
        ak_uint32 l[3]={2,3,5};
        ak_uint32 w[3]={d,d,d};
        ak_point_set(ep,ec);
        for (int cnt=1; cnt<=NUM_OF_ITER;cnt+=100,cnt-=(cnt%2==0?0:1)){
            ak_int32 kNAF[512],bases[512];
            ak_int32 i=ak_n_to_Ext_wmb_NAF(one,kNAF,bases,l,w,3,ak_mpzn512_size,ec->p,ec->r2,ec->n);
            //
            ak_mpznmax pows;
            ak_int32 A=1;
            for(int j=0;j<3;j++) A*=pow(l[j],w[j]);
            struct point Q, R_[13500],inv;
            ak_point_set_as_unit(&Q, ec );
            ak_point_set_point(R_, ep, ec );
            for (ak_int32 j=2;j<A/2;j++){
                ak_mpzn_set_ui(pows,ec->size,j);
                ak_point_pow_binary(R_+j-1,R_,pows,ec->size,ec);
            }

            //
            unsigned long long a=_rdtsc();
            for(int qwe=0;qwe<cnt;qwe++) {            for (i=i-1; i>=0; i--) {

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

                ak_point_set_point( ep, &Q, ec );
}
            unsigned long long b=_rdtsc();
            printf("Несколько оснований (основания: 2,3,5 окно:%du);%d ; %llu\n",w[0],cnt,b-a);
        }

    }
}

void Get_time(ak_curve ec,ak_mpzn512 one){
    printf("Время с учетом вычисления несмежной формы числа и предварительного набора точек\n");
    all_time(ec,one);
    all_time_mbw(ec,one);
    printf("Время с учетом вычисления предварительного набора точек\n");
    time_with_point_eval(ec,one);
    printf("Время без учета вычисления предварительного набора точек и несмежной формы\n");
    time_only_pow(ec,one);
}
#include "unistd.h"
#include "memory.h"
int main(){
    ak_mpzn512 one = ak_mpzn512_one;
    struct curve gost_e = gost_3410_2012_512_paramSetC, gost_w=gost_3410_2012_512_paramSetC_w;
    ak_curve ec = &gost_e,wc=&gost_w;
    //struct random r;
  //  ak_random_context_create_lcg(&r);
   // ak_mpzn_set_random_modulo(one, ec->q, ec->size, &r);
    ak_mpzn_set_hexstr(one,ec->size,"3D261FBD26B266A536D620C7B4BA391841B21ADE2080EDE86E95D6D76A503A58AE3DFB677357B98DD1461221FBA1E28B872C");
//    ak_random_context_destroy(&r);
    printf("Скрученные кривые Эдвардса\n");
    Get_time(ec,one);

    printf("Вейерштрасс\n");
    Get_time(wc,one);
    printf("END!\n");
    return 0;
}
