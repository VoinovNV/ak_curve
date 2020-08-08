#pragma once
#include <stdio.h>

#include <ak_mpzn.h>
#include <ak_curves.h>
typedef struct curve* ak_curve;
typedef struct point* ak_point;

struct point {
    /*! \brief x-координата точки эллиптической кривой */
    ak_uint64 x[ak_mpzn512_size];
    /*! \brief y-координата точки эллиптической кривой */
    ak_uint64 y[ak_mpzn512_size];
    /*! \brief z-координата точки эллиптической кривой */
    ak_uint64 z[ak_mpzn512_size];
};
struct ak_curve_params{
    ak_uint64 a[ak_mpzn512_size];
    ak_uint64 b[ak_mpzn512_size];

    ak_uint64 s[ak_mpzn512_size];
    ak_uint64 t[ak_mpzn512_size];
};

struct curve {
    /*! \brief � азмер параметров эллиптической кривой, исчисляемый количеством 64-х битных блоков. */
    ak_uint32 size;
    /*! \brief Кофактор эллиптической кривой - делитель порядка группы точек. */
    ak_uint32 cofactor;
    /*! \brief Коэффициент \f$ a \f$ эллиптической кривой (в представлении Монтгомери) */
    ak_uint64 a[ak_mpzn512_size];
    /*! \brief Коэффициент \f$ b \f$ эллиптической кривой (в представлении Монтгомери). */
    ak_uint64 b[ak_mpzn512_size];
    /*! \brief Модуль \f$ p \f$ эллиптической кривой. */
    ak_uint64 p[ak_mpzn512_size];
    /*! \brief Величина \f$ r^2\f$, взятая по модулю \f$ p \f$ и используемая в арифметике Монтгомери. */
    ak_uint64 r2[ak_mpzn512_size];
    /*! \brief Порядок \f$ q \f$ подгруппы, порождаемой образующей точкой \f$ P \f$. */
    ak_uint64 q[ak_mpzn512_size];
    /*! \brief Величина \f$ r^2\f$, взятая по модулю \f$ q \f$ и используемая в арифметике Монтгомери. */
    ak_uint64 r2q[ak_mpzn512_size];
    /*! \brief Точка \f$ P \f$ эллиптической кривой, порождающая подгруппу порядка \f$ q \f$. */
    struct point point;
    /*! \brief Константа \f$ n \f$, используемая в арифметике Монтгомери по модулю \f$ p \f$. */
    ak_uint64 n;
    /*! \brief Константа \f$ n_q \f$, используемая в арифметике Монтгомери по модулю \f$ q\f$. */
    ak_uint64 nq;
    /*! \brief Строка, содержащая символьную запись модуля \f$ p \f$.
     \details Используется для проверки корректного хранения параметров кривой в памяти. */
    const char *pchar;

    /*! \brief Число, указывающее форму кривой. */
    ak_int8 cform;
    /*! \brief Параметры перехода от короткой формы Вейерштрасса к форме скрученной кривой Эдвардса. */
    struct ak_curve_params edw;



};
void ak_point_set_as_unit(ak_point , ak_curve);

void ak_point_set_point(ak_point, ak_point, ak_curve);

/* ?????????? ????? ???????? ?????????? ????? ?????? */
void ak_point_set(ak_point, ak_curve);

/* ?????????? ????? ? ???????? ???? */
void ak_point_reduce(ak_point , ak_curve );

/* ???????? ????? ?? ?????????????? ??????*/
bool_t ak_point_is_ok(ak_point , ak_curve );

/* ???????? ????? ?? ?????? */
void ak_point_double(ak_point , ak_curve );

/* ???????? ???? ????? ?? ?????? */
void ak_point_add(ak_point , ak_point , ak_curve );

/* ?????????? ??????? ????? */
void ak_point_pow(ak_point , ak_point , ak_uint64 * , size_t , ak_curve );


void ak_curve_change_form(ak_curve);
/* ??????? ????? ?? ?????? ???????? ? ????? ?? ?????? ????????????*/
void ak_set_curve_form_e(ak_curve ec);
void ak_epoint_to_wpoint(ak_point wp,ak_point ep, ak_curve ec);

void ak_wpoint_to_epoint(ak_point ep, ak_point wp, ak_curve ec);

void ak_point_triple(ak_point , ak_curve );

void ak_point_quintuple(ak_point , ak_curve );

void ak_point_pow_binary( ak_point wq, ak_point wp, ak_uint64 *k, size_t size,
                          ak_curve ec );
void ak_point_pow_NAF( ak_point wq, ak_point wp, ak_int32 *kNAF,ak_int32 i,
                       ak_curve ec);
void ak_point_pow_NAF_wv(ak_point wq, ak_point wp, ak_uint64 *k, size_t size,
                         ak_curve ec);

void ak_point_pow_NAF_powof2( ak_point wq, ak_point wp,  ak_int32 *kNAF,
                             ak_int32 i,ak_uint32 w, ak_curve ec);
void ak_point_pow_NAF_powof2_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size,ak_uint32 w,
                                ak_curve ec);

void ak_point_pow_NAF_powofL( ak_point wq, ak_point wp, ak_int32 *kNAF, ak_int32 i, ak_uint32 l,ak_uint32 w,
                              ak_curve ec);
void ak_point_pow_NAF_powofL_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size, ak_uint32 l,ak_uint32 w,
                                 ak_curve ec);
void ak_point_pow_NAF_mbw( ak_point wq, ak_point wp, ak_int32 *kNAF, ak_int32 i,ak_int32 *bases, ak_uint32* l,ak_uint32* w, ak_uint8 len,
                           ak_curve ec);
void ak_point_pow_NAF_mbw_wv( ak_point wq, ak_point wp, ak_uint64 *k, size_t size, ak_uint32* l,ak_uint32* w, ak_uint8 len,
                             ak_curve ec);


ak_int32 ak_n_to_Ext_wmb_NAF(ak_uint64 *k, ak_int32* res, ak_int32* bases,
                             ak_uint32* a, ak_uint32* w, ak_int8 len,size_t size,ak_uint64 *p,ak_uint64 *r2,ak_uint64 n);

ak_int32 ak_n_to_NAF_L_w(ak_uint64 *k, ak_int32* res, ak_int32 L, ak_int32 w, size_t size, ak_uint64 *p,ak_uint64 *r2,ak_uint64 n);

ak_int32 ak_n_to_NAF_powof2(ak_uint64 *k,ak_int32 w, ak_int32* res, size_t size);
void ak_invert_p(ak_point resp,ak_point ep, ak_curve ec);
ak_int32 ak_n_to_NAF2(ak_uint64 *k, ak_int32* res, size_t size);

/* ????????? ?????????? ?????? ???????? */
const static struct curve gost_3410_2012_512_paramSetC = {
    ak_mpzn512_size, /* size */
    0,
    { 0x0000000000000239, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* a ? ????? ?????????? */
    { 0x6515a5166d05caf7, 0xae6dc7d439a723d5, 0xdc1c74edcea76671, 0x853a44eed58ae3e5, 0xc84c79f64266472e, 0xa1a4bfeccd0cf540, 0xab899e4c73783aa1, 0xde66ec2f500fc692 }, /* b ? ????? ?????????? */
	{ 0xfffffffffffffdc7, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff }, /* p */
    { 0x000000000004f0b1, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* r2 */
    { 0x94623cef47f023ed, 0xc8eda9e7a769a126, 0x4c33a9ff5147502c, 0xc98cdba46506ab00, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x3fffffffffffffff }, /* q */
    { 0xe58fa18ee6ca4eb6, 0xe79280282d956fca, 0xd016086ec2d4f903, 0x542f8f3fa490666a, 0x04f77045db49adc9, 0x314e0a57f445b20e, 0x8910352f3bea2192, 0x394c72054d8503be }, /* r2q */
    {
		{ 0x0000000000000012, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* px */
		{ 0x600303ee73001a3d, 0x905622c04b2baae7, 0xbf068c5d139732f0, 0x22dd4b650cf789ee, 0x9a56117f7b386695, 0x0fdfb0d01794368d, 0x6b99592b77a01e2a, 0x469af79d1fb1f5e1 }, /* py */
		{ 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }  /* pz */
	},
	0x58a1f7e6ce0f4c09LL, /* n */
	0x0ed9d8e0b6624e1bLL, /* nq */
    "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7",
    1,
    {
        { 0xd341ab3699869915, 0x3d6c9273ccebc4c1, 0x486b484c83cb0726, 0x9a8145b812d1a7b0, 0x2003251cadf8effa, 0x6b20d9f8b7db94f1, 0xdd0c19f57c9cc019, 0x408aa82ae77985ca }, /* a */
        { 0xb304002a3c03ce62, 0xcbe7bfdf359dc095, 0x57398fea29abadad, 0x3ce46aec38657034, 0xabf0edb5e37f775e, 0x63ccffc5280e7697, 0x6754d90e93579656, 0xc9b558b380cc6f00 }, /* b */
        { 0xa6ba96ba64be8cb4, 0x94648e0af196370a, 0x88f8e2c48c562663, 0x5eb16ec44a9d4706, 0xcdece1826f666e34, 0x9796d004ccbcc2af, 0x551d986ce321f157, 0x486644f42bfc0e5b }, /* S */
        { 0xe62e462e6780f788, 0x9d124bf8b44685f8, 0xfa04be27a2713bbd, 0x163460d278ec7b50, 0x76b769a90b110bdd, 0xf0461ffcccd77e35, 0x71ec450cbde95f1a, 0x2511275d3802a118 }, /* T */
    }

};
const static struct curve gost_3410_2012_512_paramSetC_w = {
    ak_mpzn512_size,
    4,
    { 0xd341ab3699869915, 0x3d6c9273ccebc4c1, 0x486b484c83cb0726, 0x9a8145b812d1a7b0, 0x2003251cadf8effa, 0x6b20d9f8b7db94f1, 0xdd0c19f57c9cc019, 0x408aa82ae77985ca }, /* a */
    { 0xb304002a3c03ce62, 0xcbe7bfdf359dc095, 0x57398fea29abadad, 0x3ce46aec38657034, 0xabf0edb5e37f775e, 0x63ccffc5280e7697, 0x6754d90e93579656, 0xc9b558b380cc6f00 }, /* b */
    { 0xfffffffffffffdc7, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff }, /* p */
    { 0x000000000004f0b1, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* r2 */
    { 0x94623cef47f023ed, 0xc8eda9e7a769a126, 0x4c33a9ff5147502c, 0xc98cdba46506ab00, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x3fffffffffffffff }, /* q */
    { 0xe58fa18ee6ca4eb6, 0xe79280282d956fca, 0xd016086ec2d4f903, 0x542f8f3fa490666a, 0x04f77045db49adc9, 0x314e0a57f445b20e, 0x8910352f3bea2192, 0x394c72054d8503be }, /* r2q */
    {
        { 0xc5bc7928c1950148, 0xc6fb85487eae97aa, 0xa7b9033db9ed3610, 0xa27272a7ae602bf2, 0xd385f7074cea043a, 0x2295b7a9cbaef021, 0xebe241ce593ef5de, 0xe2e31edfc23de7bd }, /* px */
        { 0xd0396e9a9addc40f, 0x04f726aa854bae07, 0xef32d85822423b63, 0xe18e2d33e3021ed2, 0x8c108c3d2090ff9b, 0x7939804d6527378b, 0xabbccff5911cb857, 0xf5ce40d95b5eb899 }, /* py */
        { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }  /* pz */
    },
    0x58a1f7e6ce0f4c09LL, /* n */
    0x0ed9d8e0b6624e1bLL, /* nq */
    "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffdc7",
    0,
    {
        { 0x0000000000000239, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }, /* a ? ????? ?????????? */
        { 0x6515a5166d05caf7, 0xae6dc7d439a723d5, 0xdc1c74edcea76671, 0x853a44eed58ae3e5, 0xc84c79f64266472e, 0xa1a4bfeccd0cf540, 0xab899e4c73783aa1, 0xde66ec2f500fc692 }, /* b ? ????? ?????????? */
        { 0xa6ba96ba64be8cb4, 0x94648e0af196370a, 0x88f8e2c48c562663, 0x5eb16ec44a9d4706, 0xcdece1826f666e34, 0x9796d004ccbcc2af, 0x551d986ce321f157, 0x486644f42bfc0e5b }, /* S */
        { 0xe62e462e6780f788, 0x9d124bf8b44685f8, 0xfa04be27a2713bbd, 0x163460d278ec7b50, 0x76b769a90b110bdd, 0xf0461ffcccd77e35, 0x71ec450cbde95f1a, 0x2511275d3802a118 }, /* T */
    }
};
