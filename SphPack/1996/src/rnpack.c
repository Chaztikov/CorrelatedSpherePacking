/* rnpack.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    real alphar, betar, cnfint, rmin, rmax, sdnrm;
} dist_;

#define dist_1 dist_

struct {
    real bound, cmptln, dtheta, eps, height;
    integer nc, npart, npart_kz__;
    real pi, porosity, rlmt;
    integer nrndsd;
    real size, theta, tolrnce;
    integer ncz;
} param_;

#define param_1 param_

struct {
    integer intensity, lshift, move;
} local_;

#define local_1 local_

struct {
    integer ncor;
    real gridln;
    integer ngrid, nzgrid, mr;
} corr_;

#define corr_1 corr_

struct {
    integer mcount, movlpcond, mpbc, mroll, nrec, ntshift, niterate;
} index_;

#define index_1 index_

struct {
    real r__[100000], rc[320000]	/* was [40][40][200] */, x[100000], 
	    xc[320000]	/* was [40][40][200] */, y[100000], yc[320000]	/* 
	    was [40][40][200] */, z__[100000], zc[320000]	/* was [40][
	    40][200] */;
    integer icnt[1600]	/* was [40][40] */;
} origdim_;

#define origdim_1 origdim_

struct {
    real base, dnot[3], finity;
} mchcom_;

#define mchcom_1 mchcom_

struct {
    integer icst2, icst3, icst4, ix, ix1, ix2, ix3, ix4, jy, jy1, jy2, jy3, 
	    jy4, kz, kz1, kz2, kz3, kz4;
    real r1, r2, r3, r4, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
} data_;

#define data_1 data_

struct {
    real rcst[12800000]	/* was [40][40][40][200] */, xcst[12800000]	/* 
	    was [40][40][40][200] */, ycst[12800000]	/* was [40][40][40][
	    200] */, zcst[12800000]	/* was [40][40][40][200] */;
    integer icntst[64000]	/* was [40][40][40] */;
} sortdim_;

#define sortdim_1 sortdim_

struct {
    integer mrolltmp;
} tmp_;

#define tmp_1 tmp_

struct {
    integer mark;
} flag_;

#define flag_1 flag_

struct {
    integer lsign;
} sign_;

#define sign_1 sign_

struct {
    real x11, x12, x13, x11max, y11, y12, y13, y11max, z11, z12, z11max, rxy, 
	    sinea, cosinea, sineb, cosineb;
} data1_;

#define data1_1 data1_

struct {
    real xshift[30], yshift[30], zshift[30];
    integer ixshift[30], jyshift[30], kzshift[30], indx[30];
} shift_;

#define shift_1 shift_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__3 = 3;
static real c_b82 = 199017.f;

/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* Main program */ int MAIN__(void)
{
    /* System generated locals */
    olist o__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void), f_open(olist *);

    /* Local variables */
    extern /* Subroutine */ int rndmpack_(void), readparam_(void), 
	    alloctcmpt_(void), sort_(void);

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 6, 0, 0, 0 };
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };


    s_wsle(&io___1);
    do_lio(&c__9, &c__1, "Start the rnpack...", (ftnlen)19);
    e_wsle();
    o__1.oerr = 0;
    o__1.ounit = 3;
    o__1.ofnmlen = 10;
    o__1.ofnm = "packin.par";
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 7;
    o__1.ofnmlen = 8;
    o__1.ofnm = "pack.par";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 8;
    o__1.ofnmlen = 9;
    o__1.ofnm = "rntmp.dat";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 9;
    o__1.ofnmlen = 8;
    o__1.ofnm = "pack.out";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_wsle(&io___2);
    do_lio(&c__9, &c__1, "Read the parameters...", (ftnlen)22);
    e_wsle();
    readparam_();
    s_wsle(&io___3);
    do_lio(&c__9, &c__1, "Generate the random packing with overlaps...", (
	    ftnlen)44);
    e_wsle();
    rndmpack_();
    s_wsle(&io___4);
    do_lio(&c__9, &c__1, "Allocate particles into compartments...", (ftnlen)
	    39);
    e_wsle();
    alloctcmpt_();
    s_wsle(&io___5);
    do_lio(&c__9, &c__1, "Seek overlap-free packing conditions...", (ftnlen)
	    39);
    e_wsle();
    sort_();
    s_wsle(&io___6);
    do_lio(&c__9, &c__1, "End the rnpack!!", (ftnlen)16);
    e_wsle();
    return 0;
} /* MAIN__ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int readparam_(void)
{
    /* Format strings */
    static char fmt_1000[] = "()";
    static char fmt_1010[] = "(/,4f9.7)";
    static char fmt_1020[] = "(/,6i7,e13.5)";
    static char fmt_1030[] = "(/,f13.7,2e13.5)";
    static char fmt_1040[] = "(/,f13.7,3i7)";

    /* Builtin functions */
    integer s_rsfe(cilist *), e_rsfe(void), do_fio(integer *, char *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 3, 0, fmt_1000, 0 };
    static cilist io___8 = { 0, 3, 0, fmt_1010, 0 };
    static cilist io___9 = { 0, 3, 0, fmt_1020, 0 };
    static cilist io___10 = { 0, 3, 0, fmt_1030, 0 };
    static cilist io___11 = { 0, 3, 0, fmt_1040, 0 };


    s_rsfe(&io___7);
    e_rsfe();
    s_rsfe(&io___8);
    do_fio(&c__1, (char *)&dist_1.alphar, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&dist_1.betar, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&dist_1.cnfint, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.porosity, (ftnlen)sizeof(real));
    e_rsfe();
    s_rsfe(&io___9);
    do_fio(&c__1, (char *)&local_1.intensity, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&index_1.mpbc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&corr_1.ncor, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&param_1.nc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&param_1.ncz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&param_1.nrndsd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&dist_1.sdnrm, (ftnlen)sizeof(real));
    e_rsfe();
    s_rsfe(&io___10);
    do_fio(&c__1, (char *)&param_1.pi, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.eps, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.tolrnce, (ftnlen)sizeof(real));
    e_rsfe();
    if (corr_1.ncor == 1) {
	s_rsfe(&io___11);
	do_fio(&c__1, (char *)&corr_1.gridln, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&corr_1.ngrid, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&corr_1.nzgrid, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&corr_1.mr, (ftnlen)sizeof(integer));
	e_rsfe();
    }
    return 0;
} /* readparam_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int rndmpack_(void)
{
    extern /* Subroutine */ int gnrtpack_(void), gnrtrdii_(void), gnrtsite_(
	    void), output_(void);

/* .....Generate the initial random packing */
    if (corr_1.ncor == 0) {
	gnrtrdii_();
	gnrtsite_();
/* .....Generate the initial random packing with spatial correlations */
    } else {
	gnrtpack_();
    }
/* .....Output the random packing with overlaps */
    output_();
    return 0;
} /* rndmpack_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int gnrtrdii_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double exp(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    extern doublereal rinvnorm_(real *);
    static integer i__;
    static real pd, vt, ravg, vpack;
    extern /* Subroutine */ int rnorm_(real *, real *, real *, real *);
    static real sample, rminlog, rmaxlog, vsphere;

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___25 = { 0, 6, 0, 0, 0 };


/* .....Set the radius range, based on the spacified confidence interval */
    r__1 = (1.f - dist_1.cnfint) * .5f;
    rminlog = dist_1.betar * rinvnorm_(&r__1) + dist_1.alphar;
    r__1 = dist_1.cnfint + (1.f - dist_1.cnfint) * .5f;
    rmaxlog = dist_1.betar * rinvnorm_(&r__1) + dist_1.alphar;
/* .....Initialize */
    param_1.npart = 0;
    vpack = 0.f;
    param_1.rlmt = .05f;
/* .....Generate the log-normal distribution of particle radii */
    dist_1.rmax = exp(rmaxlog);
    dist_1.rmin = exp(rminlog);
    param_1.cmptln = dist_1.rmax * 2.f + param_1.rlmt;
    param_1.size = param_1.cmptln * (real) param_1.nc;
    param_1.height = param_1.cmptln * (real) param_1.ncz;
    ravg = exp(dist_1.alphar);
    vt = param_1.size * param_1.size * param_1.height;
    pd = 0.f;
    if (dabs(dist_1.betar) <= param_1.tolrnce) {
	vsphere = param_1.pi * 4.f * ravg * ravg * ravg / 3.f;
	param_1.npart = (integer) (vt * (1 - param_1.porosity) / vsphere);
	i__1 = param_1.npart;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    origdim_1.r__[i__ - 1] = ravg;
	}
    } else {
	while(pd < 1 - param_1.porosity) {
	    ++param_1.npart;
	    sample = rmaxlog;
	    while(sample >= rmaxlog || sample <= rminlog) {
		rnorm_(&dist_1.alphar, &dist_1.betar, &sample, &dist_1.sdnrm);
	    }
	    origdim_1.r__[param_1.npart - 1] = exp(sample);
	    vpack += param_1.pi * 4.f * origdim_1.r__[param_1.npart - 1] * 
		    origdim_1.r__[param_1.npart - 1] * origdim_1.r__[
		    param_1.npart - 1] / 3.f;
	    pd = vpack / vt;
	}
    }
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "The maximum radius = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&dist_1.rmax, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, "The average radius = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&ravg, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___23);
    do_lio(&c__9, &c__1, "The sample size = ", (ftnlen)18);
    do_lio(&c__4, &c__1, (char *)&param_1.size, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___24);
    do_lio(&c__9, &c__1, "The sampple height = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&param_1.height, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___25);
    do_lio(&c__9, &c__1, "The number of spheres in the initial packing = ", (
	    ftnlen)47);
    do_lio(&c__3, &c__1, (char *)&param_1.npart, (ftnlen)sizeof(integer));
    e_wsle();
    param_1.bound = param_1.cmptln - dist_1.rmax * 2.f;
    return 0;
} /* gnrtrdii_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int rnorm_(real *mu, real *sigma, real *sample, real *seed)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double r_mod(real *, real *), log(doublereal), sqrt(doublereal), cos(
	    doublereal);

    /* Local variables */
    static real u1, u2, umin, umax;

    umin = 0.f;
    umax = 1.f;
L100:
    r__1 = *seed * 24298.f + 99991.f;
    u1 = r_mod(&r__1, &c_b82);
    *seed = u1;
    u1 = u1 * (umax - umin) / 199017.f + umin;
    r__1 = *seed * 24298.f + 99991.f;
    u2 = r_mod(&r__1, &c_b82);
    *seed = u2;
    u2 = u2 * (umax - umin) / 199017.f + umin;
    if (u1 == 0.f) {
	goto L100;
    }
    *sample = sqrt(log(u1) * -2.f) * cos(u2 * 6.2831799999999998f) * *sigma + 
	    *mu;
    return 0;
} /* rnorm_ */

/* ------------------------------------------------------------------------ */
doublereal rinvnorm_(real *p)
{
    /* Initialized data */

    static real connor[17] = { 8.0327350124e-17f,1.4483264644e-15f,
	    2.4668270103e-14f,3.9554295164e-13f,5.9477940136e-12f,
	    8.3507027951e-11f,1.0892221037e-9f,1.3122532964e-8f,
	    1.4503852223e-7f,1.4589169001e-6f,1.3227513228e-5f,
	    1.0683760684e-4f,7.5757575758e-4f,.0046296296296f,.02380952381f,
	    .1f,.33333333333f };
    static real rthfpi = 1.2533141373f;
    static real rrt2pi = .3989422804f;
    static real termin = 1e-11f;
    static real hstngs[6] = { 2.515517f,.802853f,.010328f,1.432788f,.189269f,
	    .001308f };

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static real a[5], b, c__, d__, e, f;
    static integer j, k, l;
    static real w, x, y, z__, p0, x0, x1, x2, x3;
    static integer ka;
    static real ord, apprxl, apprxu;


/*         INVERSE CUMULATIVE NORMAL DISTRIBUTION FUNCTION. */
/*            P = CNORM(RInvNorm(P)) */
/*     THIS IS A COPY OF ALGORITHM AS 24, */
/*     APPLIED STATISTICS, VOLUME 18, NUMBER 3, 1969 */
/*     WRITTEN BY S.W. CUNNINGHAM. */


/*         IF P IS OUT OF RANGE, RETURN POSITIVE OR */
/*     NEGATIVE INFINITY. */

    if (*p > 0.f) {
	goto L91;
    }
/* Computing 2nd power */
    r__1 = mchcom_1.base;
    ord = -(mchcom_1.finity / (r__1 * r__1));
    goto L10;
L91:
    if (*p < 1.f) {
	goto L92;
    }
    ord = mchcom_1.finity;
    goto L10;
L92:

/*         GET FIRST APPROXIMATION, X0, TO ORDINATE BY */
/*     HASTINGS' FORMULA. */

    b = *p;
    if (b > .5f) {
	b = 1.f - b;
    }
    f = -log(b);
    e = sqrt(f + f);
    x0 = -e + ((hstngs[2] * e + hstngs[1]) * e + hstngs[0]) / (((hstngs[5] * 
	    e + hstngs[4]) * e + hstngs[3]) * e + 1.f);
    if (x0 < 0.f) {
	goto L1;
    }
    x0 = 0.f;
    p0 = .5f;
    x1 = -rthfpi;
    goto L7;

/*         FIND THE AREA, P0, CORRESPONDING TO X0. */

L1:
/* Computing 2nd power */
    r__1 = x0;
    y = r__1 * r__1;
    if (x0 <= -1.9f) {
	goto L3;
    }
    y *= -.5f;

/*         (1) SERIES APPROXIMATION. */

    p0 = connor[0];
    for (l = 2; l <= 17; ++l) {
	p0 = p0 * y + connor[l - 1];
/* L2: */
    }
    p0 = (p0 * y + 1.f) * x0;
    x1 = -(p0 * rthfpi) * exp(-y);
    p0 = p0 * rrt2pi + .5f;
    goto L7;

/*         (2) CONTINUED FRACTION APPROXIMATION. */

L3:
    z__ = 1.f / y;
    a[1] = 1.f;
    a[2] = 1.f;
    a[3] = z__ + 1.f;
    a[4] = 1.f;
    w = 2.f;
L4:
    for (l = 1; l <= 3; l += 2) {
	for (j = 1; j <= 2; ++j) {
	    k = l + j;
	    ka = 7 - k;
	    a[k - 1] = a[ka - 1] + a[k - 1] * w * z__;
/* L5: */
	}
	w += 1.f;
/* L6: */
    }
    apprxu = a[1] / a[2];
    apprxl = a[4] / a[3];
    c__ = apprxu - apprxl;
    if (c__ >= termin) {
	goto L4;
    }
    x1 = apprxl / x0;
    p0 = -x1 * rrt2pi * exp(y * -.5f);

/*         GET ACCURATE VALUE OF ORDINATE BY TAYLOR SERIES. */
/*     X1, X2, AND X3 ARE DERIVATIVES FOR THE TAYLOR SERIES. */

L7:
    d__ = f + log(p0);
    x2 = x0 * x1 * x1 - x1;
/* Computing 3rd power */
    r__1 = x1;
    x3 = r__1 * (r__1 * r__1) + x0 * 2.f * x1 * x2 - x2;
    x = ((x3 * d__ / 3.f + x2) * d__ / 2.f + x1) * d__ + x0;
    ord = x;
    if (*p > .5f) {
	ord = -x;
    }
L10:
    ret_val = ord;
    return ret_val;
} /* rinvnorm_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int gnrtsite_(void)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;
    extern doublereal ran_(integer *);

    i__1 = param_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	origdim_1.x[i__ - 1] = (r__1 = param_1.size * ran_(&param_1.nrndsd), 
		dabs(r__1));
	origdim_1.y[i__ - 1] = (r__1 = param_1.size * ran_(&param_1.nrndsd), 
		dabs(r__1));
	origdim_1.z__[i__ - 1] = (r__1 = param_1.height * ran_(&
		param_1.nrndsd), dabs(r__1));
    }
    return 0;
} /* gnrtsite_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int gnrtpack_(void)
{
    /* Format strings */
    static char fmt_2000[] = "(2(2x,e13.5))";
    static char fmt_1000[] = "(///)";
    static char fmt_1010[] = "(10x,f10.4)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    olist o__1;

    /* Builtin functions */
    integer f_open(olist *);
    double exp(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsfe(cilist *), e_rsfe(void), s_wsle(cilist *), do_lio(integer 
	    *, integer *, char *, ftnlen), e_wsle(void), s_rsle(cilist *), 
	    e_rsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer ncradius[100];
    extern doublereal rinvnorm_(real *);
    static integer k, l;
    static real pd, pf;
    static integer ix, jy, kz;
    static real vt;
    extern doublereal ran_(integer *);
    static real cpdf[100], delr, ravg, pfun[1000000]	/* was [100][100][100]
	     */, vpack;
    extern /* Subroutine */ int rnorm_(real *, real *, real *, real *);
    static real sample, radius;
    static integer lsignal, isample, nradius[100];
    static real rminlog, rmaxlog, xsample, ysample, zsample;

    /* Fortran I/O blocks */
    static cilist io___72 = { 0, 12, 0, fmt_2000, 0 };
    static cilist io___73 = { 0, 4, 0, fmt_1000, 0 };
    static cilist io___77 = { 0, 4, 0, fmt_1010, 0 };
    static cilist io___84 = { 0, 6, 0, 0, 0 };
    static cilist io___85 = { 0, 6, 0, 0, 0 };
    static cilist io___86 = { 0, 6, 0, 0, 0 };
    static cilist io___87 = { 0, 6, 0, 0, 0 };
    static cilist io___88 = { 0, 6, 0, 0, 0 };
    static cilist io___89 = { 0, 6, 0, 0, 0 };
    static cilist io___90 = { 0, 6, 0, 0, 0 };
    static cilist io___91 = { 0, 5, 0, 0, 0 };


    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = 10;
    o__1.ofnm = "pfield.dat";
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 12;
    o__1.ofnmlen = 8;
    o__1.ofnm = "cpsd.out";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/* .....Generate the cumulative radius size distribution */
    r__1 = (1.f - dist_1.cnfint) * .5f;
    rminlog = dist_1.betar * rinvnorm_(&r__1) + dist_1.alphar;
    r__1 = dist_1.cnfint + (1.f - dist_1.cnfint) * .5f;
    rmaxlog = dist_1.betar * rinvnorm_(&r__1) + dist_1.alphar;
    dist_1.rmax = exp(rmaxlog);
    dist_1.rmin = exp(rminlog);
    ravg = exp(dist_1.alphar);
    delr = (dist_1.rmax - dist_1.rmin) / (real) corr_1.mr;
    param_1.rlmt = .05f;
    param_1.cmptln = dist_1.rmax * 2.f + param_1.rlmt;
    param_1.size = param_1.cmptln * (real) param_1.nc;
    param_1.height = param_1.cmptln * (real) param_1.ncz;
/* .....Initialize */
    i__1 = corr_1.mr;
    for (k = 1; k <= i__1; ++k) {
	nradius[k - 1] = 0;
    }
    param_1.npart = 0;
    pd = 0.f;
    vpack = 0.f;
    vt = param_1.size * param_1.size * param_1.height;
    for (isample = 1; isample <= 100000; ++isample) {
	sample = rmaxlog;
	while(sample >= rmaxlog || sample <= rminlog) {
	    rnorm_(&dist_1.alphar, &dist_1.betar, &sample, &dist_1.sdnrm);
	}
	radius = exp(sample);
	i__1 = corr_1.mr;
	for (k = 1; k <= i__1; ++k) {
	    if (radius >= dist_1.rmin + delr * (k - 1) && radius < 
		    dist_1.rmin + delr * k) {
		++nradius[k - 1];
	    }
	}
    }
    ncradius[0] = nradius[0];
    i__1 = corr_1.mr - 1;
    for (k = 1; k <= i__1; ++k) {
	ncradius[k] = nradius[k] + ncradius[k - 1];
    }
/* .....Normalize the distribution */
    i__1 = corr_1.mr;
    for (k = 1; k <= i__1; ++k) {
	cpdf[k - 1] = (real) ncradius[k - 1] / (real) ncradius[corr_1.mr - 1];
	s_wsfe(&io___72);
	r__1 = dist_1.rmin + delr * (k - 1);
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cpdf[k - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
/* .....Read in the random P-field data */
    s_rsfe(&io___73);
    e_rsfe();
    i__1 = corr_1.nzgrid;
    for (kz = 1; kz <= i__1; ++kz) {
	i__2 = corr_1.ngrid;
	for (jy = 1; jy <= i__2; ++jy) {
	    i__3 = corr_1.ngrid;
	    for (ix = 1; ix <= i__3; ++ix) {
		s_rsfe(&io___77);
		do_fio(&c__1, (char *)&pfun[ix + (jy + kz * 100) * 100 - 
			10101], (ftnlen)sizeof(real));
		e_rsfe();
/*    				write(*,*) ix,jy,kz,Pfun(ix,jy,kz) */
	    }
	}
    }
/* .....Create the initial random packing with spatial correlations */
    while(pd < 1 - param_1.porosity) {
	++param_1.npart;
	xsample = (r__1 = param_1.size * ran_(&param_1.nrndsd), dabs(r__1));
	i__1 = corr_1.ngrid;
	for (l = 1; l <= i__1; ++l) {
	    if (xsample >= corr_1.gridln * (l - 1) && xsample < corr_1.gridln 
		    * l) {
		ix = l;
	    }
	}
	ysample = (r__1 = param_1.size * ran_(&param_1.nrndsd), dabs(r__1));
	i__1 = corr_1.ngrid;
	for (l = 1; l <= i__1; ++l) {
	    if (ysample >= corr_1.gridln * (l - 1) && ysample < corr_1.gridln 
		    * l) {
		jy = l;
	    }
	}
	zsample = (r__1 = param_1.height * ran_(&param_1.nrndsd), dabs(r__1));
	i__1 = corr_1.nzgrid;
	for (l = 1; l <= i__1; ++l) {
	    if (zsample >= corr_1.gridln * (l - 1) && zsample < corr_1.gridln 
		    * l) {
		kz = l;
	    }
	}
	origdim_1.x[param_1.npart - 1] = xsample;
	origdim_1.y[param_1.npart - 1] = ysample;
	origdim_1.z__[param_1.npart - 1] = zsample;
	pf = pfun[ix + (jy + kz * 100) * 100 - 10101];
	if (pf < cpdf[0]) {
	    origdim_1.r__[param_1.npart - 1] = dist_1.rmin - delr * ran_(&
		    param_1.nrndsd);
	} else if (pf >= cpdf[corr_1.mr - 1]) {
	    origdim_1.r__[param_1.npart - 1] = dist_1.rmax + delr * ran_(&
		    param_1.nrndsd);
	} else {
	    i__1 = corr_1.mr - 1;
	    for (k = 1; k <= i__1; ++k) {
		if (pf >= cpdf[k - 1] && pf < cpdf[k]) {
		    origdim_1.r__[param_1.npart - 1] = dist_1.rmin + delr * (
			    k - 1) + delr * ran_(&param_1.nrndsd);
		}
	    }
	}
	vpack += param_1.pi * 4.f * origdim_1.r__[param_1.npart - 1] * 
		origdim_1.r__[param_1.npart - 1] * origdim_1.r__[
		param_1.npart - 1] / 3.f;
	pd = vpack / vt;
    }
    s_wsle(&io___84);
    do_lio(&c__9, &c__1, "The maximum radius = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&dist_1.rmax, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___85);
    do_lio(&c__9, &c__1, "The minimum radius = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&dist_1.rmin, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___86);
    do_lio(&c__9, &c__1, "The average radius = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&ravg, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___87);
    do_lio(&c__9, &c__1, "The sample size = ", (ftnlen)18);
    do_lio(&c__4, &c__1, (char *)&param_1.size, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___88);
    do_lio(&c__9, &c__1, "The sampple height = ", (ftnlen)21);
    do_lio(&c__4, &c__1, (char *)&param_1.height, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___89);
    do_lio(&c__9, &c__1, "The number of spheres in the initial packing = ", (
	    ftnlen)47);
    do_lio(&c__3, &c__1, (char *)&param_1.npart, (ftnlen)sizeof(integer));
    e_wsle();
    param_1.bound = param_1.cmptln - dist_1.rmax * 2.f;
    s_wsle(&io___90);
    do_lio(&c__9, &c__1, "Do you want to start rearranging process(1--Yes, 0"
	    "--No)?", (ftnlen)56);
    e_wsle();
    s_rsle(&io___91);
    do_lio(&c__3, &c__1, (char *)&lsignal, (ftnlen)sizeof(integer));
    e_rsle();
    if (lsignal == 0) {
	s_stop("", (ftnlen)0);
    }
    return 0;
} /* gnrtpack_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int output_(void)
{
    /* Format strings */
    static char fmt_2010[] = "(4(2x,e13.5))";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    static cilist io___94 = { 0, 8, 0, fmt_2010, 0 };


/*     WRITE(8,2000) */
/* 		WRITE(8,*) (x(i), y(i), z(i), r(i), i = 1, npart) */
    i__1 = param_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_wsfe(&io___94);
	do_fio(&c__1, (char *)&origdim_1.x[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.y[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.z__[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.r__[i__ - 1], (ftnlen)sizeof(real));
	e_wsfe();
    }
/* L2000: */
    return 0;
} /* output_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int alloctcmpt_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(4(2x,e13.5))";
    static char fmt_2000[] = "(\002npart_kz\002,/,\002ix  jy  ncpart\002,/"
	    ",\002xc  yc  zc  rc\002)";
    static char fmt_2010[] = "(i7)";
    static char fmt_2020[] = "(3i7)";
    static char fmt_2030[] = "(4(2x,e13.5))";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    alist al__1;

    /* Builtin functions */
    integer f_rew(alist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(void), s_wsfe(cilist *), e_wsfe(void);

    /* Local variables */
    static integer i__, ic, mz;
    extern /* Subroutine */ int allocating_(integer *, integer *), 
	    initialcnt_(void);
    static integer ncpart;

    /* Fortran I/O blocks */
    static cilist io___96 = { 0, 8, 0, fmt_1000, 0 };
    static cilist io___97 = { 0, 9, 0, fmt_2000, 0 };
    static cilist io___99 = { 0, 9, 0, fmt_2010, 0 };
    static cilist io___101 = { 0, 9, 0, fmt_2020, 0 };
    static cilist io___103 = { 0, 9, 0, fmt_2030, 0 };


/* .....Read in particles and allocate them into compartments */
    al__1.aerr = 0;
    al__1.aunit = 8;
    f_rew(&al__1);
/* 		READ(8,*) */
/* 		READ(8,*) (x(i), y(i), z(i), r(i), i = 1, npart) */
    i__1 = param_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_rsfe(&io___96);
	do_fio(&c__1, (char *)&origdim_1.x[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.y[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.z__[i__ - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&origdim_1.r__[i__ - 1], (ftnlen)sizeof(real));
	e_rsfe();
    }
    s_wsfe(&io___97);
    e_wsfe();
    mz = 0;
    i__1 = param_1.ncz;
    for (data_1.kz = 1; data_1.kz <= i__1; ++data_1.kz) {
	initialcnt_();
	i__2 = param_1.npart;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    allocating_(&i__, &mz);
	}
/* .....Sum up the total counts in the layer */
	i__2 = param_1.nc;
	for (data_1.jy = 1; data_1.jy <= i__2; ++data_1.jy) {
	    i__3 = param_1.nc;
	    for (data_1.ix = 1; data_1.ix <= i__3; ++data_1.ix) {
		param_1.npart_kz__ += origdim_1.icnt[data_1.ix + data_1.jy * 
			40 - 41];
	    }
	}
/* .....Write particles in file after allocation */
	s_wsfe(&io___99);
	do_fio(&c__1, (char *)&param_1.npart_kz__, (ftnlen)sizeof(integer));
	e_wsfe();
	i__2 = param_1.nc;
	for (data_1.jy = 1; data_1.jy <= i__2; ++data_1.jy) {
	    i__3 = param_1.nc;
	    for (data_1.ix = 1; data_1.ix <= i__3; ++data_1.ix) {
		ncpart = origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41];
		s_wsfe(&io___101);
		do_fio(&c__1, (char *)&data_1.ix, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&data_1.jy, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ncpart, (ftnlen)sizeof(integer));
		e_wsfe();
/* 					WRITE(9,*) (xc(ix,jy,ic), yc(ix,jy,ic), zc(ix,jy,ic), */
/*     1                     rc(ix,jy,ic), ic = 1, ncpart) */
		i__4 = ncpart;
		for (ic = 1; ic <= i__4; ++ic) {
		    s_wsfe(&io___103);
		    do_fio(&c__1, (char *)&origdim_1.xc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.yc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.zc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.rc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    e_wsfe();
		}
	    }
	}
	++mz;
    }
    return 0;
} /* alloctcmpt_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int initialcnt_(void)
{
    /* System generated locals */
    integer i__1, i__2;

/* .....Initiallize counts in each compartment */
    i__1 = param_1.nc;
    for (data_1.jy = 1; data_1.jy <= i__1; ++data_1.jy) {
	i__2 = param_1.nc;
	for (data_1.ix = 1; data_1.ix <= i__2; ++data_1.ix) {
	    origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41] = 0;
	}
    }
/* .....Initialize count in each layer */
    param_1.npart_kz__ = 0;
    return 0;
} /* initialcnt_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int allocating_(integer *i__, integer *mz)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer ic, my, mx;

    my = 0;
    i__1 = param_1.nc;
    for (data_1.jy = 1; data_1.jy <= i__1; ++data_1.jy) {
	mx = 0;
	i__2 = param_1.nc;
	for (data_1.ix = 1; data_1.ix <= i__2; ++data_1.ix) {
	    if (origdim_1.x[*i__ - 1] >= mx * param_1.cmptln && origdim_1.x[*
		    i__ - 1] < (mx + 1) * param_1.cmptln && origdim_1.y[*i__ 
		    - 1] >= my * param_1.cmptln && origdim_1.y[*i__ - 1] < (
		    my + 1) * param_1.cmptln && origdim_1.z__[*i__ - 1] >= *
		    mz * param_1.cmptln && origdim_1.z__[*i__ - 1] < (*mz + 1)
		     * param_1.cmptln) {
		++origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41];
		ic = origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41];
		origdim_1.xc[data_1.ix + (data_1.jy + ic * 40) * 40 - 1641] = 
			origdim_1.x[*i__ - 1];
		origdim_1.yc[data_1.ix + (data_1.jy + ic * 40) * 40 - 1641] = 
			origdim_1.y[*i__ - 1];
		origdim_1.zc[data_1.ix + (data_1.jy + ic * 40) * 40 - 1641] = 
			origdim_1.z__[*i__ - 1];
		origdim_1.rc[data_1.ix + (data_1.jy + ic * 40) * 40 - 1641] = 
			origdim_1.r__[*i__ - 1];
	    }
	    ++mx;
	}
	++my;
    }
    return 0;
} /* allocating_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int sort_(void)
{
    /* Format strings */
    static char fmt_1000[] = "(//)";
    static char fmt_1010[] = "(4(2x,e13.5))";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    alist al__1;

    /* Builtin functions */
    integer f_rew(alist *), s_rsfe(cilist *), e_rsfe(void), s_rsle(cilist *), 
	    do_lio(integer *, integer *, char *, ftnlen), e_rsle(void), 
	    do_fio(integer *, char *, ftnlen), s_wsle(cilist *), e_wsle(void);

    /* Local variables */
    extern /* Subroutine */ int finalset_(void), findzmin_(integer *, integer 
	    *);
    static integer ic, lz, ic1;
    extern /* Subroutine */ int initialsort_(integer *);
    static integer ipart;
    static real dummy;
    static integer ncpart, ncpart1;
    extern /* Subroutine */ int reorder_(integer *, integer *), sorting_(
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___108 = { 0, 9, 0, fmt_1000, 0 };
    static cilist io___109 = { 0, 9, 0, 0, 0 };
    static cilist io___110 = { 0, 9, 0, 0, 0 };
    static cilist io___114 = { 0, 9, 0, fmt_1010, 0 };
    static cilist io___115 = { 0, 6, 0, 0, 0 };


/* .....Initialize after-sort counts */
    i__1 = param_1.ncz + 1;
    for (lz = 1; lz <= i__1; ++lz) {
	initialsort_(&lz);
    }
    index_1.niterate = 0;
/* .....Sort particles layer by layer */
    al__1.aerr = 0;
    al__1.aunit = 9;
    f_rew(&al__1);
    s_rsfe(&io___108);
    e_rsfe();
    i__1 = param_1.ncz;
    for (lz = 1; lz <= i__1; ++lz) {
	s_rsle(&io___109);
	do_lio(&c__3, &c__1, (char *)&param_1.npart_kz__, (ftnlen)sizeof(
		integer));
	e_rsle();
	i__2 = param_1.nc;
	for (data_1.jy = 1; data_1.jy <= i__2; ++data_1.jy) {
	    i__3 = param_1.nc;
	    for (data_1.ix = 1; data_1.ix <= i__3; ++data_1.ix) {
		s_rsle(&io___110);
		do_lio(&c__4, &c__1, (char *)&dummy, (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&dummy, (ftnlen)sizeof(real));
		do_lio(&c__3, &c__1, (char *)&ncpart, (ftnlen)sizeof(integer))
			;
		e_rsle();
		origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41] = ncpart;
		i__4 = ncpart;
		for (ic = 1; ic <= i__4; ++ic) {
		    s_rsfe(&io___114);
		    do_fio(&c__1, (char *)&origdim_1.xc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.yc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.zc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    do_fio(&c__1, (char *)&origdim_1.rc[data_1.ix + (
			    data_1.jy + ic * 40) * 40 - 1641], (ftnlen)sizeof(
			    real));
		    e_rsfe();
		}
	    }
	}
	s_wsle(&io___115);
	do_lio(&c__9, &c__1, "Working on the vertical cell: ", (ftnlen)30);
	do_lio(&c__3, &c__1, (char *)&lz, (ftnlen)sizeof(integer));
	e_wsle();
	i__2 = param_1.npart_kz__;
	for (ipart = 1; ipart <= i__2; ++ipart) {
	    findzmin_(&ic1, &ncpart1);
	    reorder_(&ic1, &ncpart1);
	    sorting_(&lz);
	}
    }
/* .....Calculate the final porosity and write out the packing data */
    finalset_();
    return 0;
} /* sort_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int initialsort_(integer *lz)
{
    /* System generated locals */
    integer i__1, i__2;

/* .....Initialize counts in a layer of lz after sort */
    i__1 = param_1.nc;
    for (data_1.jy = 1; data_1.jy <= i__1; ++data_1.jy) {
	i__2 = param_1.nc;
	for (data_1.ix = 1; data_1.ix <= i__2; ++data_1.ix) {
	    sortdim_1.icntst[data_1.ix + (data_1.jy + *lz * 40) * 40 - 1641] =
		     0;
	}
    }
    return 0;
} /* initialsort_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findzmin_(integer *ic1, integer *ncpart1)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer ic, ncpart;

    data_1.z1 = param_1.size;
    i__1 = param_1.nc;
    for (data_1.jy = 1; data_1.jy <= i__1; ++data_1.jy) {
	i__2 = param_1.nc;
	for (data_1.ix = 1; data_1.ix <= i__2; ++data_1.ix) {
	    ncpart = origdim_1.icnt[data_1.ix + data_1.jy * 40 - 41];
	    i__3 = ncpart;
	    for (ic = 1; ic <= i__3; ++ic) {
		if (origdim_1.zc[data_1.ix + (data_1.jy + ic * 40) * 40 - 
			1641] < data_1.z1) {
		    data_1.z1 = origdim_1.zc[data_1.ix + (data_1.jy + ic * 40)
			     * 40 - 1641];
		    *ic1 = ic;
		    *ncpart1 = ncpart;
		    data_1.ix1 = data_1.ix;
		    data_1.jy1 = data_1.jy;
		}
	    }
	}
    }
    return 0;
} /* findzmin_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int reorder_(integer *ic1, integer *ncpart1)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ic;

    data_1.x1 = origdim_1.xc[data_1.ix1 + (data_1.jy1 + *ic1 * 40) * 40 - 
	    1641];
    data_1.y1 = origdim_1.yc[data_1.ix1 + (data_1.jy1 + *ic1 * 40) * 40 - 
	    1641];
    data_1.z1 = origdim_1.zc[data_1.ix1 + (data_1.jy1 + *ic1 * 40) * 40 - 
	    1641];
    data_1.r1 = origdim_1.rc[data_1.ix1 + (data_1.jy1 + *ic1 * 40) * 40 - 
	    1641];
    --origdim_1.icnt[data_1.ix1 + data_1.jy1 * 40 - 41];
    i__1 = *ncpart1;
    for (ic = *ic1 + 1; ic <= i__1; ++ic) {
	origdim_1.xc[data_1.ix1 + (data_1.jy1 + (ic - 1) * 40) * 40 - 1641] = 
		origdim_1.xc[data_1.ix1 + (data_1.jy1 + ic * 40) * 40 - 1641];
	origdim_1.yc[data_1.ix1 + (data_1.jy1 + (ic - 1) * 40) * 40 - 1641] = 
		origdim_1.yc[data_1.ix1 + (data_1.jy1 + ic * 40) * 40 - 1641];
	origdim_1.zc[data_1.ix1 + (data_1.jy1 + (ic - 1) * 40) * 40 - 1641] = 
		origdim_1.zc[data_1.ix1 + (data_1.jy1 + ic * 40) * 40 - 1641];
	origdim_1.rc[data_1.ix1 + (data_1.jy1 + (ic - 1) * 40) * 40 - 1641] = 
		origdim_1.rc[data_1.ix1 + (data_1.jy1 + ic * 40) * 40 - 1641];
    }
    return 0;
} /* reorder_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int sorting_(integer *lz)
{
    extern /* Subroutine */ int register_(void), findlocbc_(integer *), 
	    findloc_(integer *), findxyz_(void);

    data_1.kz1 = *lz;
    index_1.ntshift = 0;
/* .....Check overlap conditions */
    index_1.mroll = 0;
L10:
    findxyz_();
    local_1.lshift = 1;
    if (data_1.z1 >= param_1.bound && data_1.x1 >= param_1.bound && data_1.y1 
	    >= param_1.bound && data_1.x1 <= param_1.size - param_1.bound && 
	    data_1.y1 <= param_1.size - param_1.bound) {
	findloc_(lz);
	if (index_1.nrec == 1) {
	    goto L20;
	}
    }
    if (data_1.z1 >= param_1.bound && (data_1.x1 < param_1.bound || data_1.x1 
	    > param_1.size - param_1.bound || data_1.y1 < param_1.bound || 
	    data_1.y1 > param_1.size - param_1.bound)) {
	findlocbc_(lz);
	if (index_1.mcount == 11) {
	    goto L10;
	}
    }
L20:
    register_();
    return 0;
} /* sorting_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findxyz_(void)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static integer ixx, jyy, kzz, mflag;
    static real zorig;
    extern /* Subroutine */ int shift1_(integer *, integer *, integer *, 
	    integer *), shift2_(integer *, integer *, integer *, integer *), 
	    shift3_(integer *, integer *, integer *, integer *), shift4_(
	    integer *, integer *, integer *, integer *);

    mflag = 1;
    index_1.movlpcond = 0;
    tmp_1.mrolltmp = 0;
    zorig = data_1.z1;
    while(mflag == 1) {
	if (data_1.ix1 == param_1.nc) {
	    ixx = param_1.nc - 1;
	}
	if (data_1.ix1 < param_1.nc) {
	    ixx = data_1.ix1;
	}
	if (data_1.jy1 == param_1.nc) {
	    jyy = param_1.nc - 1;
	}
	if (data_1.jy1 < param_1.nc) {
	    jyy = data_1.jy1;
	}
	if (data_1.kz1 == param_1.ncz + 1) {
	    kzz = param_1.ncz;
	}
	if (data_1.kz1 < param_1.ncz + 1) {
	    kzz = data_1.kz1;
	}
	mflag = 0;
	if (data_1.kz1 == 1 || data_1.kz1 == param_1.ncz + 1) {
	    if (data_1.ix1 == 1 && data_1.jy1 == 1 || data_1.ix1 == 
		    param_1.nc && data_1.jy1 == 1 || data_1.ix1 == 1 && 
		    data_1.jy1 == param_1.nc || data_1.ix1 == param_1.nc && 
		    data_1.jy1 == param_1.nc) {
		shift1_(&ixx, &jyy, &kzz, &mflag);
	    } else if (data_1.ix1 == 1 || data_1.ix1 == param_1.nc || 
		    data_1.jy1 == 1 || data_1.jy1 == param_1.nc) {
		shift2_(&ixx, &jyy, &kzz, &mflag);
	    } else {
		shift3_(&ixx, &jyy, &kzz, &mflag);
	    }
	} else if (data_1.jy1 == 1 || data_1.jy1 == param_1.nc) {
	    if (data_1.ix1 == 1 || data_1.ix1 == param_1.nc) {
		shift2_(&ixx, &jyy, &kzz, &mflag);
	    } else {
		shift3_(&ixx, &jyy, &kzz, &mflag);
	    }
	} else if (data_1.ix1 == 1 || data_1.ix1 == param_1.nc) {
	    shift3_(&ixx, &jyy, &kzz, &mflag);
	} else {
	    shift4_(&ixx, &jyy, &kzz, &mflag);
	}
    }
    if (index_1.mroll == 0 && (r__1 = zorig - data_1.z1, dabs(r__1)) > 
	    param_1.tolrnce) {
	index_1.mroll = 1;
    }
    if (tmp_1.mrolltmp == 3 && index_1.movlpcond == 0) {
	index_1.mroll = tmp_1.mrolltmp;
    }
    return 0;
} /* findxyz_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int shift1_(integer *ixx, integer *jyy, integer *kzz, 
	integer *mflag)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int chkovlp_(integer *);

    i__1 = *kzz + 1;
    for (data_1.kz = *kzz; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *jyy + 1;
	for (data_1.jy = *jyy; data_1.jy <= i__2; ++data_1.jy) {
	    i__3 = *ixx + 1;
	    for (data_1.ix = *ixx; data_1.ix <= i__3; ++data_1.ix) {
		chkovlp_(mflag);
		if (*mflag == 1 || index_1.movlpcond == 1) {
		    goto L10;
		}
	    }
	}
    }
L10:
    return 0;
} /* shift1_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int shift2_(integer *ixx, integer *jyy, integer *kzz, 
	integer *mflag)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int chkovlp_(integer *);

    if (data_1.ix1 >= 2 && data_1.ix1 < param_1.nc) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx - 1; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    } else if (data_1.jy1 >= 2 && data_1.jy1 < param_1.nc) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy - 1; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    } else if (data_1.kz1 >= 2 && data_1.kz1 <= param_1.ncz) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz - 1; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    }
L10:
    return 0;
} /* shift2_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int shift3_(integer *ixx, integer *jyy, integer *kzz, 
	integer *mflag)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int chkovlp_(integer *);

    if (data_1.ix1 == 1 || data_1.ix1 == param_1.nc) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz - 1; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy - 1; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    } else if (data_1.jy1 == 1 || data_1.jy1 == param_1.nc) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz - 1; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx - 1; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    } else if (data_1.kz1 == 1 || data_1.kz1 == param_1.ncz + 1) {
	i__1 = *kzz + 1;
	for (data_1.kz = *kzz; data_1.kz <= i__1; ++data_1.kz) {
	    i__2 = *jyy + 1;
	    for (data_1.jy = *jyy - 1; data_1.jy <= i__2; ++data_1.jy) {
		i__3 = *ixx + 1;
		for (data_1.ix = *ixx - 1; data_1.ix <= i__3; ++data_1.ix) {
		    chkovlp_(mflag);
		    if (*mflag == 1 || index_1.movlpcond == 1) {
			goto L10;
		    }
		}
	    }
	}
    }
L10:
    return 0;
} /* shift3_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int shift4_(integer *ixx, integer *jyy, integer *kzz, 
	integer *mflag)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int chkovlp_(integer *);

    i__1 = *kzz + 1;
    for (data_1.kz = *kzz - 1; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *jyy + 1;
	for (data_1.jy = *jyy - 1; data_1.jy <= i__2; ++data_1.jy) {
	    i__3 = *ixx + 1;
	    for (data_1.ix = *ixx - 1; data_1.ix <= i__3; ++data_1.ix) {
		chkovlp_(mflag);
		if (*mflag == 1 || index_1.movlpcond == 1) {
		    goto L10;
		}
	    }
	}
    }
L10:
    return 0;
} /* shift4_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int chkovlp_(integer *mflag)
{
    extern /* Subroutine */ int roll1loc_(void), roll2loc_(void), findzmax_(
	    integer *);

    if (index_1.mroll == 0) {
	findzmax_(mflag);
    } else if (index_1.mroll == 1) {
	roll1loc_();
	goto L10;
    } else if (index_1.mroll == 2) {
	roll2loc_();
    }
L10:
    return 0;
} /* chkovlp_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findzmax_(integer *mflag)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real d__, rr, xr, yr, zr, zp, dxy;
    static integer icst;
    static real sumr, sinea;
    static integer ncntst;
    static real cosinea;

    ncntst = sortdim_1.icntst[data_1.ix + (data_1.jy + data_1.kz * 40) * 40 - 
	    1641];
    if (ncntst == 0) {
	goto L20;
    }
    i__1 = ncntst;
    for (icst = 1; icst <= i__1; ++icst) {
	xr = sortdim_1.xcst[data_1.ix + (data_1.jy + (data_1.kz + icst * 40) *
		 40) * 40 - 65641];
	yr = sortdim_1.ycst[data_1.ix + (data_1.jy + (data_1.kz + icst * 40) *
		 40) * 40 - 65641];
	zr = sortdim_1.zcst[data_1.ix + (data_1.jy + (data_1.kz + icst * 40) *
		 40) * 40 - 65641];
	rr = sortdim_1.rcst[data_1.ix + (data_1.jy + (data_1.kz + icst * 40) *
		 40) * 40 - 65641];
	d__ = sqrt((data_1.x1 - xr) * (data_1.x1 - xr) + (data_1.y1 - yr) * (
		data_1.y1 - yr) + (data_1.z1 - zr) * (data_1.z1 - zr));
	sumr = data_1.r1 + rr;
	if (d__ < sumr - param_1.eps) {
L10:
	    dxy = sqrt((data_1.x1 - xr) * (data_1.x1 - xr) + (data_1.y1 - yr) 
		    * (data_1.y1 - yr));
	    if (dxy <= param_1.tolrnce) {
		data_1.x1 += 1e-4f;
		data_1.y1 += 1e-4f;
		if (data_1.x1 >= param_1.cmptln * data_1.ix1 && data_1.ix1 < 
			param_1.nc) {
		    ++data_1.ix1;
		}
		if (data_1.y1 >= param_1.cmptln * data_1.jy1 && data_1.jy1 < 
			param_1.nc) {
		    ++data_1.jy1;
		}
		goto L10;
	    }
	    cosinea = dxy / sumr;
	    if (cosinea > 1.f) {
		cosinea = 1.f;
	    }
	    sinea = sqrt(1.f - cosinea * cosinea);
	    zp = sumr * sinea;
	    data_1.z1 = zr + zp;
	    data_1.x2 = xr;
	    data_1.y2 = yr;
	    data_1.z2 = zr;
	    data_1.r2 = rr;
	    data_1.ix2 = data_1.ix;
	    data_1.jy2 = data_1.jy;
	    data_1.kz2 = data_1.kz;
	    data_1.icst2 = icst;
	    if (data_1.z1 >= param_1.cmptln * data_1.kz && data_1.kz1 <= 
		    param_1.ncz) {
		data_1.kz1 = data_1.kz + 1;
	    }
	    *mflag = 1;
	    goto L20;
	}
    }
L20:
    return 0;
} /* findzmax_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int roll1loc_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icst, ncntst;
    extern /* Subroutine */ int calr1loc_(integer *);

    ncntst = sortdim_1.icntst[data_1.ix + (data_1.jy + data_1.kz * 40) * 40 - 
	    1641];
    i__1 = ncntst;
    for (icst = 1; icst <= i__1; ++icst) {
	if (! (data_1.ix == data_1.ix2 && data_1.jy == data_1.jy2 && 
		data_1.kz == data_1.kz2 && icst == data_1.icst2)) {
	    calr1loc_(&icst);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* roll1loc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int calr1loc_(integer *icst)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real d__, rr, xr, yr, zr, sumr;

    xr = sortdim_1.xcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    yr = sortdim_1.ycst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    zr = sortdim_1.zcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    rr = sortdim_1.rcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    d__ = sqrt((data_1.x1 - xr) * (data_1.x1 - xr) + (data_1.y1 - yr) * (
	    data_1.y1 - yr) + (data_1.z1 - zr) * (data_1.z1 - zr));
    sumr = data_1.r1 + rr;
    if (d__ >= sumr - param_1.eps && d__ <= sumr) {
	data_1.x3 = xr;
	data_1.y3 = yr;
	data_1.z3 = zr;
	data_1.r3 = rr;
	data_1.ix3 = data_1.ix;
	data_1.jy3 = data_1.jy;
	data_1.kz3 = data_1.kz;
	data_1.icst3 = *icst;
	index_1.mroll = 2;
    } else if (d__ < sumr - param_1.eps) {
	if (param_1.dtheta <= param_1.eps) {
	    data_1.x3 = xr;
	    data_1.y3 = yr;
	    data_1.z3 = zr;
	    data_1.r3 = rr;
	    data_1.ix3 = data_1.ix;
	    data_1.jy3 = data_1.jy;
	    data_1.kz3 = data_1.kz;
	    data_1.icst3 = *icst;
	    index_1.mroll = 2;
	} else {
	    param_1.theta -= param_1.dtheta;
	    param_1.dtheta *= .5f;
	    index_1.movlpcond = 1;
	}
    }
    return 0;
} /* calr1loc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int roll2loc_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icst, ncntst;
    extern /* Subroutine */ int calr2loc_(integer *);

    ncntst = sortdim_1.icntst[data_1.ix + (data_1.jy + data_1.kz * 40) * 40 - 
	    1641];
    i__1 = ncntst;
    for (icst = 1; icst <= i__1; ++icst) {
	if (! (data_1.ix == data_1.ix2 && data_1.jy == data_1.jy2 && 
		data_1.kz == data_1.kz2 && icst == data_1.icst2) && ! (
		data_1.ix == data_1.ix3 && data_1.jy == data_1.jy3 && 
		data_1.kz == data_1.kz3 && icst == data_1.icst3)) {
	    calr2loc_(&icst);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* roll2loc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int calr2loc_(integer *icst)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real d__, rr, xr, yr, zr, sumr;

    xr = sortdim_1.xcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    yr = sortdim_1.ycst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    zr = sortdim_1.zcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    rr = sortdim_1.rcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    d__ = sqrt((data_1.x1 - xr) * (data_1.x1 - xr) + (data_1.y1 - yr) * (
	    data_1.y1 - yr) + (data_1.z1 - zr) * (data_1.z1 - zr));
    sumr = data_1.r1 + rr;
    if (d__ >= sumr - param_1.eps && d__ <= sumr) {
	data_1.x4 = xr;
	data_1.y4 = yr;
	data_1.z4 = zr;
	data_1.r4 = rr;
	data_1.ix4 = data_1.ix;
	data_1.jy4 = data_1.jy;
	data_1.kz4 = data_1.kz;
	data_1.icst4 = *icst;
	tmp_1.mrolltmp = 3;
    } else if (d__ < sumr - param_1.eps) {
	if (param_1.dtheta <= param_1.eps) {
	    data_1.x4 = xr;
	    data_1.y4 = yr;
	    data_1.z4 = zr;
	    data_1.r4 = rr;
	    data_1.ix4 = data_1.ix;
	    data_1.jy4 = data_1.jy;
	    data_1.kz4 = data_1.kz;
	    data_1.icst4 = *icst;
	    tmp_1.mrolltmp = 3;
	} else {
	    param_1.theta -= param_1.dtheta;
	    param_1.dtheta *= .5f;
	    index_1.movlpcond = 1;
	}
    }
    return 0;
} /* calr2loc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int register_(void)
{
    static integer icst1;

    icst1 = sortdim_1.icntst[data_1.ix1 + (data_1.jy1 + data_1.kz1 * 40) * 40 
	    - 1641] + 1;
    sortdim_1.icntst[data_1.ix1 + (data_1.jy1 + data_1.kz1 * 40) * 40 - 1641] 
	    = icst1;
    sortdim_1.xcst[data_1.ix1 + (data_1.jy1 + (data_1.kz1 + icst1 * 40) * 40) 
	    * 40 - 65641] = data_1.x1;
    sortdim_1.ycst[data_1.ix1 + (data_1.jy1 + (data_1.kz1 + icst1 * 40) * 40) 
	    * 40 - 65641] = data_1.y1;
    sortdim_1.zcst[data_1.ix1 + (data_1.jy1 + (data_1.kz1 + icst1 * 40) * 40) 
	    * 40 - 65641] = data_1.z1;
    sortdim_1.rcst[data_1.ix1 + (data_1.jy1 + (data_1.kz1 + icst1 * 40) * 40) 
	    * 40 - 65641] = data_1.r1;
/*      write(*,2010) xcst(ix1,jy1,kz1,icst1), ycst(ix1,jy1,kz1,icst1), */
/*     1              zcst(ix1,jy1,kz1,icst1), rcst(ix1,jy1,kz1,icst1) */
/*      write(*,2020) ix1,jy1,kz1,icst1 */
/*      write(*,2000) rcst(ix1,jy1,kz1,icst1), ix1,jy1,kz1,icst1 */
/* L2000: */
/* L2010: */
/* L2020: */
    return 0;
} /* register_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int shiftdown_(void)
{
    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static real delz, ravg;
    static integer ncdiv;
    extern /* Subroutine */ int findxyz_(void);

    ravg = exp(dist_1.alphar);
    ncdiv = 30;
    delz = ravg / (real) ncdiv;
    while(data_1.z1 >= param_1.bound) {
	data_1.z1 -= delz;
	if (data_1.z1 < param_1.cmptln * (data_1.kz1 - 1) && data_1.kz1 > 1) {
	    --data_1.kz1;
	}
	index_1.mroll = 0;
	findxyz_();
	if (index_1.mroll == 1) {
	    goto L10;
	}
    }
L10:
    return 0;
} /* shiftdown_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findlocbc_(integer *lz)
{
    extern /* Subroutine */ int rollitbc_(void), shiftdown_(void), 
	    relocatebc_(void);

    index_1.mcount = 0;
    index_1.nrec = 0;
    index_1.mroll = 0;
    while(index_1.nrec != 1) {
	if (index_1.mroll == 0) {
	    shiftdown_();
	}
	if (data_1.z1 < param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	}
	if (index_1.mroll == 1) {
	    if ((data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		    param_1.bound) && (data_1.y1 < param_1.bound || data_1.y1 
		    > param_1.size - param_1.bound)) {
		if (data_1.kz1 > *lz) {
		    goto L10;
		} else {
		    index_1.nrec = 1;
		    goto L20;
		}
	    } else {
		rollitbc_();
	    }
	}
	if (index_1.mroll >= 2) {
	    ++index_1.niterate;
	    if (data_1.kz1 > *lz) {
		goto L10;
	    } else {
		index_1.nrec = 1;
		goto L20;
	    }
	}
L10:
	if (data_1.kz1 > *lz) {
	    relocatebc_();
	}
	if (index_1.mcount == 11) {
	    goto L20;
	}
    }
L20:
    return 0;
} /* findlocbc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int rollitbc_(void)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal), sin(doublereal), cos(
	    doublereal);

    /* Local variables */
    extern /* Subroutine */ int chkbcmpt_(real *, integer *);
    static real o1, o2, p1, p2;
    extern /* Subroutine */ int findxyzbc_(real *, integer *);
    static real d12, o12;
    static integer lp1;
    static real delp, p1tmp, z1tmp, angle, sinet;
    static integer lp1tmp, kz1tmp;
    static real rtrack, dptort, cosinet;

    if (data_1.x1 < param_1.bound || data_1.x1 > param_1.size - param_1.bound)
	     {
	flag_1.mark = 1;
	o1 = data_1.x1;
	o2 = data_1.x2;
	p1 = data_1.y1;
	p2 = data_1.y2;
	lp1 = data_1.jy1;
    } else {
	flag_1.mark = 2;
	o1 = data_1.y1;
	o2 = data_1.y2;
	p1 = data_1.x1;
	p2 = data_1.x2;
	lp1 = data_1.ix1;
    }
    if ((r__1 = p1 - p2, dabs(r__1)) <= param_1.tolrnce) {
	p1 += 1e-4f;
	if (p1 >= param_1.cmptln * lp1 && lp1 < param_1.nc) {
	    ++lp1;
	}
	rtrack = data_1.z1 - data_1.z2;
	sinet = 1e-4f / rtrack;
	if (sinet > 1.f) {
	    sinet = 1.f;
	}
	cosinet = sqrt(1.f - sinet * sinet);
	data_1.z1 = data_1.z2 + rtrack * cosinet;
    }
    d12 = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2) + (data_1.z1 - data_1.z2) *
	     (data_1.z1 - data_1.z2));
    o12 = (r__1 = o1 - o2, dabs(r__1));
    rtrack = sqrt(d12 * d12 - o12 * o12);
    delp = (r__1 = p1 - p2, dabs(r__1));
    dptort = delp / rtrack;
    if (dptort > 1.f) {
	param_1.theta = param_1.pi * .5f;
    } else {
	param_1.theta = asin(dptort);
    }
    param_1.dtheta = .15f;
    angle = param_1.pi * .5f;
    while(param_1.theta < angle) {
	p1tmp = p1;
	z1tmp = data_1.z1;
	lp1tmp = lp1;
	kz1tmp = data_1.kz1;
L10:
	param_1.theta += param_1.dtheta;
	if (p1 > p2) {
	    p1 = p2 + rtrack * sin(param_1.theta);
	} else {
	    p1 = p2 - rtrack * sin(param_1.theta);
	}
	data_1.z1 = data_1.z2 + rtrack * cos(param_1.theta);
	chkbcmpt_(&p1, &lp1);
	findxyzbc_(&p1, &lp1);
	if (index_1.movlpcond == 1) {
	    goto L10;
	}
	if (index_1.mroll == 2 || data_1.z1 < param_1.bound) {
	    if (flag_1.mark == 1) {
		if (param_1.dtheta <= param_1.eps) {
		    data_1.y1 = p1tmp;
		    data_1.z1 = z1tmp;
		    data_1.jy1 = lp1tmp;
		    data_1.kz1 = kz1tmp;
		} else {
		    data_1.y1 = p1;
		    data_1.jy1 = lp1;
		}
	    } else {
		if (param_1.dtheta <= param_1.eps) {
		    data_1.x1 = p1tmp;
		    data_1.z1 = z1tmp;
		    data_1.ix1 = lp1tmp;
		    data_1.kz1 = kz1tmp;
		} else {
		    data_1.x1 = p1;
		    data_1.ix1 = lp1;
		}
	    }
	    index_1.nrec = 1;
	    goto L20;
	}
    }
    index_1.mroll = 0;
    if (flag_1.mark == 1) {
	data_1.y1 = p1;
	data_1.jy1 = lp1;
    } else {
	data_1.x1 = p1;
	data_1.ix1 = lp1;
    }
L20:
    return 0;
} /* rollitbc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int chkbcmpt_(real *p1, integer *lp1)
{
    if (*p1 >= param_1.cmptln * *lp1 && *lp1 < param_1.nc) {
	++(*lp1);
    }
    if (*p1 < param_1.cmptln * (*lp1 - 1) && *lp1 >= 2) {
	--(*lp1);
    }
    if (data_1.z1 < param_1.cmptln * (data_1.kz1 - 1) && data_1.kz1 >= 2) {
	--data_1.kz1;
    }
    return 0;
} /* chkbcmpt_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findxyzbc_(real *p1, integer *lp1)
{
    static integer lpp, kzp;
    extern /* Subroutine */ int bcroll1_(real *, integer *, integer *, 
	    integer *), bcroll2_(real *, integer *, integer *, integer *), 
	    bcroll3_(real *, integer *, integer *, integer *), bcroll4_(real *
	    , integer *, integer *, integer *);

    index_1.movlpcond = 0;
    if (*lp1 == param_1.nc) {
	lpp = param_1.nc - 1;
    }
    if (*lp1 < param_1.nc) {
	lpp = *lp1;
    }
    if (data_1.kz1 == param_1.ncz + 1) {
	kzp = param_1.ncz;
    }
    if (data_1.kz1 < param_1.ncz + 1) {
	kzp = data_1.kz1;
    }
    if (*lp1 == 1 && data_1.kz1 == 1 || *lp1 == param_1.nc && data_1.kz1 == 1 
	    || *lp1 == 1 && data_1.kz1 == param_1.ncz + 1 || *lp1 == 
	    param_1.nc && data_1.kz1 == param_1.ncz + 1) {
	bcroll1_(p1, lp1, &lpp, &kzp);
    } else if (*lp1 == 1 || *lp1 == param_1.nc) {
	bcroll2_(p1, lp1, &lpp, &kzp);
    } else if (data_1.kz1 == 1 || data_1.kz1 == param_1.ncz + 1) {
	bcroll3_(p1, lp1, &lpp, &kzp);
    } else {
	bcroll4_(p1, lp1, &lpp, &kzp);
    }
    return 0;
} /* findxyzbc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int bcroll1_(real *p1, integer *lp1, integer *lpp, integer *
	kzp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int chkbcovlp_(real *, integer *, integer *);
    static integer lp;

    i__1 = *kzp + 1;
    for (data_1.kz = *kzp; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *lpp + 1;
	for (lp = *lpp; lp <= i__2; ++lp) {
	    chkbcovlp_(p1, lp1, &lp);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* bcroll1_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int bcroll2_(real *p1, integer *lp1, integer *lpp, integer *
	kzp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int chkbcovlp_(real *, integer *, integer *);
    static integer lp;

    i__1 = *kzp + 1;
    for (data_1.kz = *kzp - 1; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *lpp + 1;
	for (lp = *lpp; lp <= i__2; ++lp) {
	    chkbcovlp_(p1, lp1, &lp);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* bcroll2_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int bcroll3_(real *p1, integer *lp1, integer *lpp, integer *
	kzp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int chkbcovlp_(real *, integer *, integer *);
    static integer lp;

    i__1 = *kzp + 1;
    for (data_1.kz = *kzp; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *lpp + 1;
	for (lp = *lpp - 1; lp <= i__2; ++lp) {
	    chkbcovlp_(p1, lp1, &lp);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* bcroll3_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int bcroll4_(real *p1, integer *lp1, integer *lpp, integer *
	kzp)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int chkbcovlp_(real *, integer *, integer *);
    static integer lp;

    i__1 = *kzp + 1;
    for (data_1.kz = *kzp - 1; data_1.kz <= i__1; ++data_1.kz) {
	i__2 = *lpp + 1;
	for (lp = *lpp - 1; lp <= i__2; ++lp) {
	    chkbcovlp_(p1, lp1, &lp);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* bcroll4_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int chkbcovlp_(real *p1, integer *lp1, integer *lp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int calrlocbc_(integer *);
    static integer icst, ncntst;

    if (flag_1.mark == 1) {
	data_1.y1 = *p1;
	data_1.jy = *lp;
	data_1.ix = data_1.ix1;
    } else {
	data_1.x1 = *p1;
	data_1.ix = *lp;
	data_1.jy = data_1.jy1;
    }
    ncntst = sortdim_1.icntst[data_1.ix + (data_1.jy + data_1.kz * 40) * 40 - 
	    1641];
    i__1 = ncntst;
    for (icst = 1; icst <= i__1; ++icst) {
	if (! (data_1.ix == data_1.ix2 && data_1.jy == data_1.jy2 && 
		data_1.kz == data_1.kz2 && icst == data_1.icst2)) {
	    calrlocbc_(&icst);
	    if (index_1.movlpcond == 1) {
		goto L10;
	    }
	}
    }
L10:
    return 0;
} /* chkbcovlp_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int calrlocbc_(integer *icst)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real d__, rr, xr, yr, zr, sumr;

    xr = sortdim_1.xcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    yr = sortdim_1.ycst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    zr = sortdim_1.zcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    rr = sortdim_1.rcst[data_1.ix + (data_1.jy + (data_1.kz + *icst * 40) * 
	    40) * 40 - 65641];
    d__ = sqrt((data_1.x1 - xr) * (data_1.x1 - xr) + (data_1.y1 - yr) * (
	    data_1.y1 - yr) + (data_1.z1 - zr) * (data_1.z1 - zr));
    sumr = data_1.r1 + rr;
    if (d__ >= sumr - param_1.eps && d__ <= sumr) {
	index_1.mroll = 2;
    } else if (d__ < sumr - param_1.eps) {
	if (param_1.dtheta <= param_1.eps) {
	    index_1.mroll = 2;
	} else {
	    param_1.theta -= param_1.dtheta;
	    param_1.dtheta *= .5f;
	    index_1.movlpcond = 1;
	}
    }
    return 0;
} /* calrlocbc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int relocatebc_(void)
{
    static integer nchalf;

    nchalf = (integer) (param_1.nc * .5f);
    if (index_1.mcount == 10) {
	if (flag_1.mark == 1) {
	    if (data_1.x1 < param_1.bound) {
		data_1.x1 += param_1.cmptln;
		if (data_1.x1 >= param_1.cmptln) {
		    ++data_1.ix1;
		}
	    } else {
		data_1.x1 -= param_1.cmptln;
		if (data_1.x1 < param_1.size - param_1.cmptln) {
		    --data_1.ix1;
		}
	    }
	} else {
	    if (data_1.y1 < param_1.bound) {
		data_1.y1 += param_1.cmptln;
		if (data_1.y1 >= param_1.cmptln) {
		    ++data_1.jy1;
		}
	    } else {
		data_1.y1 -= param_1.cmptln;
		if (data_1.y1 < param_1.size - param_1.cmptln) {
		    --data_1.jy1;
		}
	    }
	}
	goto L20;
    }
    if (index_1.ntshift >= 60) {
	data_1.r1 *= .9f;
	index_1.mcount = 0;
	index_1.ntshift = 0;
    }
L10:
    if (index_1.mcount == 0) {
	if ((data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		param_1.bound) && (data_1.y1 < param_1.bound || data_1.y1 > 
		param_1.size - param_1.bound)) {
	    if (data_1.jy1 <= nchalf) {
		sign_1.lsign = 1;
	    } else {
		sign_1.lsign = -1;
	    }
	} else if (data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		param_1.bound) {
	    if (data_1.jy1 <= nchalf) {
		sign_1.lsign = 1;
	    } else {
		sign_1.lsign = -1;
	    }
	} else {
	    if (data_1.ix1 <= nchalf) {
		sign_1.lsign = 1;
	    } else {
		sign_1.lsign = -1;
	    }
	}
    }
    if (flag_1.mark == 1) {
	if (data_1.y1 > param_1.size) {
	    data_1.y1 = param_1.size * 2.f - data_1.y1;
	} else if (data_1.y1 < 0.f) {
	    data_1.y1 = -data_1.y1;
	}
	data_1.y1 += param_1.cmptln * sign_1.lsign;
	data_1.jy1 += sign_1.lsign;
    } else {
	if (data_1.x1 > param_1.size) {
	    data_1.x1 = param_1.size * 2.f - data_1.x1;
	} else if (data_1.x1 < 0.f) {
	    data_1.x1 = -data_1.x1;
	}
	data_1.x1 += param_1.cmptln * sign_1.lsign;
	data_1.ix1 += sign_1.lsign;
    }
    if (data_1.ix1 == 0 || data_1.ix1 == param_1.nc + 1) {
	data_1.x1 -= param_1.cmptln * 2.f * (real) sign_1.lsign;
	data_1.ix1 -= sign_1.lsign << 1;
	data_1.r1 *= .9f;
	index_1.mcount = 0;
	index_1.ntshift = 0;
	goto L10;
    } else if (data_1.jy1 == 0 || data_1.jy1 == param_1.nc + 1) {
	data_1.y1 -= param_1.cmptln * 2.f * (real) sign_1.lsign;
	data_1.jy1 -= sign_1.lsign << 1;
	data_1.r1 *= .9f;
	index_1.mcount = 0;
	index_1.ntshift = 0;
	goto L10;
    }
L20:
    ++index_1.mcount;
    index_1.mroll = 0;
    index_1.nrec = 0;
    ++index_1.ntshift;
    return 0;
} /* relocatebc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int findloc_(integer *lz)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    extern /* Subroutine */ int relocate_(void), shiftdown_(void), 
	    settlerlmt_(void);
    static integer nloop;
    extern /* Subroutine */ int tmprec_(void), rollon1_(void), rollon2_(void),
	     chkcond_(void);

    index_1.mcount = 0;
    index_1.nrec = 0;
L10:
    if ((r__1 = data_1.r1 - param_1.rlmt, dabs(r__1)) <= param_1.tolrnce) {
	while(index_1.nrec != 1) {
	    if (index_1.mroll == 0) {
		shiftdown_();
	    }
	    if (data_1.z1 < param_1.bound) {
		index_1.nrec = 1;
		goto L20;
	    }
	    if (index_1.nrec == 2) {
		goto L20;
	    }
	    if (index_1.mroll == 2) {
		settlerlmt_();
		++index_1.niterate;
	    }
	}
    } else {
	nloop = 0;
	while(index_1.nrec != 1) {
	    if (index_1.mroll == 0) {
		shiftdown_();
	    }
	    if (data_1.z1 < param_1.bound) {
		index_1.nrec = 1;
		goto L20;
	    }
	    if (index_1.mroll == 1) {
		rollon1_();
	    }
	    if (index_1.nrec == 2) {
		goto L20;
	    }
	    if (index_1.mroll == 2) {
		rollon2_();
	    }
	    if (index_1.mroll == 3) {
		chkcond_();
	    }
	    if (index_1.nrec != 1) {
		++nloop;
	    }
	    if (nloop >= 30) {
		relocate_();
		nloop = 0;
		local_1.lshift = 1;
	    }
	    if (index_1.nrec == 1) {
		++index_1.niterate;
		if (data_1.z1 >= param_1.bound) {
		    tmprec_();
		}
	    }
	    if (index_1.nrec == 1 && local_1.lshift == local_1.move && 
		    data_1.kz1 > *lz) {
		relocate_();
		nloop = 0;
		local_1.lshift = 1;
	    }
	    if ((r__1 = data_1.r1 - param_1.rlmt, dabs(r__1)) <= 
		    param_1.tolrnce) {
		goto L10;
	    }
	}
    }
L20:
    return 0;
} /* findloc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int rollon1_(void)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal), sin(doublereal), cos(
	    doublereal);

    /* Local variables */
    static real d12, dx, dxy, beta, angle;
    extern /* Subroutine */ int chkcmpt_(void);
    static real dxytod12;
    extern /* Subroutine */ int findxyz_(void);
    static real dxtodxy;

    dx = (r__1 = data_1.x1 - data_1.x2, dabs(r__1));
    dxy = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2));
    dxtodxy = dx / dxy;
    if (dxtodxy > 1.f) {
	beta = param_1.pi * .5f;
    } else {
	beta = asin(dxtodxy);
    }
    d12 = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2) + (data_1.z1 - data_1.z2) *
	     (data_1.z1 - data_1.z2));
    dxytod12 = dxy / d12;
    if (dxytod12 > 1.f) {
	dxytod12 = 1.f;
    }
    if (data_1.z1 > data_1.z2) {
	param_1.theta = asin(dxytod12);
    } else {
	index_1.mroll = 0;
	goto L20;
    }
    param_1.dtheta = .15f;
    angle = param_1.pi * .5f;
    while(param_1.theta < angle) {
L10:
	param_1.theta += param_1.dtheta;
	dxy = d12 * sin(param_1.theta);
	if (data_1.x1 >= data_1.x2) {
	    data_1.x1 = data_1.x2 + dxy * sin(beta);
	} else {
	    data_1.x1 = data_1.x2 - dxy * sin(beta);
	}
	if (data_1.y1 >= data_1.y2) {
	    data_1.y1 = data_1.y2 + dxy * cos(beta);
	} else {
	    data_1.y1 = data_1.y2 - dxy * cos(beta);
	}
	data_1.z1 = data_1.z2 + d12 * cos(param_1.theta);
	chkcmpt_();
	index_1.mroll = 1;
	findxyz_();
	if (index_1.movlpcond == 1) {
	    goto L10;
	}
	if (data_1.z1 < param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	} else if (data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		param_1.bound || data_1.y1 < param_1.bound || data_1.y1 > 
		param_1.size - param_1.bound) {
	    if (data_1.z1 <= data_1.z2) {
		index_1.mroll = 0;
	    }
	    index_1.nrec = 2;
	    goto L20;
	} else if (index_1.mroll >= 2) {
	    goto L20;
	}
    }
    index_1.mroll = 0;
L20:
    return 0;
} /* rollon1_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int chkcmpt_(void)
{
    if (data_1.x1 >= param_1.cmptln * data_1.ix1 && data_1.ix1 < param_1.nc) {
	++data_1.ix1;
    }
    if (data_1.x1 < param_1.cmptln * (data_1.ix1 - 1) && data_1.ix1 >= 2) {
	--data_1.ix1;
    }
    if (data_1.y1 >= param_1.cmptln * data_1.jy1 && data_1.jy1 < param_1.nc) {
	++data_1.jy1;
    }
    if (data_1.y1 < param_1.cmptln * (data_1.jy1 - 1) && data_1.jy1 >= 2) {
	--data_1.jy1;
    }
    if (data_1.z1 >= param_1.cmptln * data_1.kz1 && data_1.kz1 <= param_1.ncz)
	     {
	++data_1.kz1;
    }
    if (data_1.z1 < param_1.cmptln * (data_1.kz1 - 1) && data_1.kz1 >= 2) {
	--data_1.kz1;
    }
    return 0;
} /* chkcmpt_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int settlerlmt_(void)
{
    extern /* Subroutine */ int calrlmtloc_(real *, real *, real *);
    static real alpha, angle, alpha2, alpha3;
    extern /* Subroutine */ int chkcmpt_(void), findxyz_(void);

    calrlmtloc_(&alpha, &alpha2, &alpha3);
    chkcmpt_();
    index_1.mroll = 0;
    findxyz_();
    angle = param_1.pi * .5f;
    if (index_1.mroll == 0 && alpha + alpha3 <= angle) {
	index_1.nrec = 1;
    } else if (index_1.mroll == 0 && alpha + alpha3 > angle) {
	data_1.x2 = data_1.x3;
	data_1.y2 = data_1.y3;
	data_1.z2 = data_1.z3;
	data_1.r2 = data_1.r3;
	data_1.ix2 = data_1.ix3;
	data_1.jy2 = data_1.jy3;
	data_1.kz2 = data_1.kz3;
	data_1.icst2 = data_1.icst3;
	index_1.mroll = 1;
    }
    return 0;
} /* settlerlmt_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int calrlmtloc_(real *alpha, real *alpha2, real *alpha3)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), asin(doublereal), acos(doublereal), sin(
	    doublereal), cos(doublereal);

    /* Local variables */
    static real d12, d13, d23, dx1, dy1, dz1, ang, dz12, dz13, ry23, beta, 
	    dxy12, dxy13, dxy23;
    static integer nsign;
    static real ratio2, ratio3, z23tod23;

    d23 = sqrt((data_1.x2 - data_1.x3) * (data_1.x2 - data_1.x3) + (data_1.y2 
	    - data_1.y3) * (data_1.y2 - data_1.y3) + (data_1.z2 - data_1.z3) *
	     (data_1.z2 - data_1.z3));
    d12 = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2) + (data_1.z1 - data_1.z2) *
	     (data_1.z1 - data_1.z2));
    d13 = sqrt((data_1.x1 - data_1.x3) * (data_1.x1 - data_1.x3) + (data_1.y1 
	    - data_1.y3) * (data_1.y1 - data_1.y3) + (data_1.z1 - data_1.z3) *
	     (data_1.z1 - data_1.z3));
    dxy23 = sqrt(d23 * d23 - (data_1.z2 - data_1.z3) * (data_1.z2 - data_1.z3)
	    );
    z23tod23 = (r__1 = data_1.z2 - data_1.z3, dabs(r__1)) / d23;
    if (z23tod23 > 1.f) {
	z23tod23 = 1.f;
    }
    *alpha = asin(z23tod23);
    ratio2 = (d12 * d12 + d23 * d23 - d13 * d13) * .5f / (d12 * d23);
    ratio3 = (d13 * d13 + d23 * d23 - d12 * d12) * .5f / (d13 * d23);
    if (ratio2 > 1.f) {
	ratio2 = 1.f;
    }
    if (ratio3 > 1.f) {
	ratio3 = 1.f;
    }
    *alpha2 = acos(ratio2);
    *alpha3 = acos(ratio3);
    ry23 = (r__1 = data_1.y2 - data_1.y3, dabs(r__1)) / dxy23;
    if ((r__1 = 1.f - ry23, dabs(r__1)) <= param_1.eps) {
	beta = param_1.pi * .5f;
    } else {
	beta = asin(ry23);
    }
    ang = param_1.pi * .5f;
    if (data_1.z2 >= data_1.z3) {
	dz1 = d12 * sin(*alpha2 - *alpha);
	data_1.z1 = data_1.z2 + dz1;
	dz13 = data_1.z1 - data_1.z3;
	dxy13 = sqrt(d13 * d13 - dz13 * dz13);
	dx1 = dxy13 * cos(beta);
	dy1 = dxy13 * sin(beta);
    } else {
	dz1 = d13 * sin(*alpha3 - *alpha);
	data_1.z1 = data_1.z3 + dz1;
	dz12 = data_1.z1 - data_1.z2;
	dxy12 = sqrt(d12 * d12 - dz12 * dz12);
	if (*alpha + *alpha2 <= ang) {
	    dx1 = (dxy23 - dxy12) * cos(beta);
	    dy1 = (dxy23 - dxy12) * sin(beta);
	} else {
	    dx1 = (dxy23 + dxy12) * cos(beta);
	    dy1 = (dxy23 + dxy12) * sin(beta);
	}
    }
    if (data_1.z2 >= data_1.z3) {
	if (*alpha + *alpha3 > ang) {
	    nsign = -1;
	} else {
	    nsign = 1;
	}
    } else {
	nsign = 1;
    }
    if (data_1.x2 >= data_1.x3 && data_1.y2 >= data_1.y3) {
	data_1.x1 = data_1.x3 + nsign * dx1;
	data_1.y1 = data_1.y3 + nsign * dy1;
    }
    if (data_1.x2 >= data_1.x3 && data_1.y2 < data_1.y3) {
	data_1.x1 = data_1.x3 + nsign * dx1;
	data_1.y1 = data_1.y3 - nsign * dy1;
    }
    if (data_1.x2 < data_1.x3 && data_1.y2 >= data_1.y3) {
	data_1.x1 = data_1.x3 - nsign * dx1;
	data_1.y1 = data_1.y3 + nsign * dy1;
    }
    if (data_1.x2 < data_1.x3 && data_1.y2 < data_1.y3) {
	data_1.x1 = data_1.x3 - nsign * dx1;
	data_1.y1 = data_1.y3 - nsign * dy1;
    }
    return 0;
} /* calrlmtloc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int rollon2_(void)
{
    /* System generated locals */
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    extern /* Subroutine */ int vertroll_(real *, real *, real *);
    static real x0, y0, z0, d12, d13, rm, rn;
    extern /* Subroutine */ int calrlmtloc_(real *, real *, real *);
    static real xtmp, ytmp, ztmp, x1max, y1max, z1max, alpha, angle, ratio, 
	    alpha2, alpha3, rtrack;
    extern /* Subroutine */ int calrloc_(real *, real *, real *, real *, real 
	    *, real *, real *);

/* .....Find the highest point of the rolling track by using subroutine */
/*     CalRlmtLoc */
    xtmp = data_1.x1;
    ytmp = data_1.y1;
    ztmp = data_1.z1;
    calrlmtloc_(&alpha, &alpha2, &alpha3);
    x1max = data_1.x1;
    y1max = data_1.y1;
    z1max = data_1.z1;
    data_1.x1 = xtmp;
    data_1.y1 = ytmp;
    data_1.z1 = ztmp;
/* .....Find the center and radius of the track */
    d13 = sqrt((data_1.x1 - data_1.x3) * (data_1.x1 - data_1.x3) + (data_1.y1 
	    - data_1.y3) * (data_1.y1 - data_1.y3) + (data_1.z1 - data_1.z3) *
	     (data_1.z1 - data_1.z3));
    d12 = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2) + (data_1.z1 - data_1.z2) *
	     (data_1.z1 - data_1.z2));
    rtrack = d13 * sin(alpha3);
    rn = d13 * cos(alpha3);
    rm = d12 * cos(alpha2);
    ratio = rm / rn;
    angle = param_1.pi * .5f;
    if ((r__1 = alpha3 - angle, dabs(r__1)) <= param_1.tolrnce) {
	x0 = data_1.x3;
	y0 = data_1.y3;
	z0 = data_1.z3;
    } else {
	x0 = (data_1.x2 + ratio * data_1.x3) / (ratio + 1);
	y0 = (data_1.y2 + ratio * data_1.y3) / (ratio + 1);
	z0 = (data_1.z2 + ratio * data_1.z3) / (ratio + 1);
    }
/* .....Start rolling */
    if ((r__1 = data_1.z2 - data_1.z3, dabs(r__1)) <= .01f) {
/* 		IF(ABS(z2-z3).LE.eps) THEN */
	vertroll_(&x1max, &y1max, &z1max);
    } else {
	calrloc_(&x1max, &y1max, &z1max, &x0, &y0, &z0, &rtrack);
    }
    return 0;
} /* rollon2_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int vertroll_(real *x1max, real *y1max, real *z1max)
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double asin(doublereal), cos(doublereal), sin(doublereal), atan(
	    doublereal), sqrt(doublereal);

    /* Local variables */
    static real dx, dxy, beta, beta1, angle, sinet, rtrack;
    extern /* Subroutine */ int chkcmpt_(void), findxyz_(void);
    static real dxtodxy, dxytort;

    rtrack = *z1max - data_1.z2;
    if ((r__1 = *z1max - data_1.z1, dabs(r__1)) <= param_1.tolrnce) {
	if ((r__1 = data_1.y2 - data_1.y3, dabs(r__1)) <= param_1.tolrnce) {
	    data_1.x1 += 1e-5f;
	    if (data_1.x1 >= param_1.cmptln * data_1.ix1 && data_1.ix1 < 
		    param_1.nc) {
		++data_1.ix1;
	    }
	    beta = 0.f;
	    sinet = 1e-5f / rtrack;
	    if (sinet > 1.f) {
		sinet = 1.f;
	    }
	    param_1.theta = asin(sinet);
	    data_1.z1 = data_1.z2 + rtrack * cos(param_1.theta);
	} else {
	    param_1.theta = 1e-6f;
	    dxy = rtrack * sin(param_1.theta);
	    beta1 = atan((r__1 = data_1.x3 - data_1.x2, dabs(r__1)) / (r__2 = 
		    data_1.y3 - data_1.y2, dabs(r__2)));
	    beta = param_1.pi * .5f - beta1;
	    dx = dxy * sin(beta);
	}
    } else {
	dx = (r__1 = data_1.x1 - *x1max, dabs(r__1));
	dxy = sqrt((data_1.x1 - *x1max) * (data_1.x1 - *x1max) + (data_1.y1 - 
		*y1max) * (data_1.y1 - *y1max));
	dxtodxy = dx / dxy;
	if (dxtodxy > 1.f) {
	    dxtodxy = 1.f;
	}
	beta = asin(dxtodxy);
	dxytort = dxy / rtrack;
	if (dxytort > 1.f) {
	    dxytort = 1.f;
	}
	param_1.theta = asin(dxy / rtrack);
    }
    param_1.dtheta = .15f;
    angle = param_1.pi * .5f;
    while(param_1.theta < angle) {
L10:
	param_1.theta += param_1.dtheta;
	dxy = rtrack * sin(param_1.theta);
	if (data_1.x1 >= *x1max) {
	    data_1.x1 = *x1max + dxy * sin(beta);
	} else {
	    data_1.x1 = *x1max - dxy * sin(beta);
	}
	if (data_1.y1 >= *y1max) {
	    data_1.y1 = *y1max + dxy * cos(beta);
	} else {
	    data_1.y1 = *y1max - dxy * cos(beta);
	}
	data_1.z1 = data_1.z2 + rtrack * cos(param_1.theta);
	chkcmpt_();
	index_1.mroll = 2;
	findxyz_();
	if (index_1.movlpcond == 1) {
	    goto L10;
	}
	if (data_1.z1 < param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	} else if (data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		param_1.bound || data_1.y1 < param_1.bound || data_1.y1 > 
		param_1.size - param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	} else if (index_1.mroll == 3) {
	    goto L20;
	}
    }
    index_1.mroll = 0;
L20:
    return 0;
} /* vertroll_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int calrloc_(real *x1max, real *y1max, real *z1max, real *x0,
	 real *y0, real *z0, real *rtrack)
{
    /* Builtin functions */
    double acos(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    extern /* Subroutine */ int transbwd_(real *, real *, real *), transfwd_(
	    real *, real *, real *, real *, real *, real *, real *);
    static real x1tmp, y1tmp, z1tmp, angle;
    static integer msign, ix1tmp, jy1tmp, kz1tmp;
    static real x13tort;
    extern /* Subroutine */ int chkcmpt_(void), findxyz_(void);

/* .....Transform axes */
    transfwd_(x1max, y1max, z1max, x0, y0, z0, rtrack);
/* .....Roll and check overlaps */
    x13tort = data1_1.x13 / *rtrack;
    if (x13tort > 1.f) {
	x13tort = 1.f;
    }
    if (x13tort < -1.f) {
	x13tort = -1.f;
    }
    param_1.theta = acos(data1_1.x13 / *rtrack);
    if (data1_1.y13 >= 0.f) {
	msign = 1;
    } else {
	msign = -1;
    }
    param_1.dtheta = .15f;
    angle = param_1.pi * .5f;
    if (param_1.theta > angle) {
	index_1.mroll = 0;
	goto L20;
    }
    while(param_1.theta < angle) {
	x1tmp = data_1.x1;
	y1tmp = data_1.y1;
	z1tmp = data_1.z1;
	ix1tmp = data_1.ix1;
	jy1tmp = data_1.jy1;
	kz1tmp = data_1.kz1;
L10:
	param_1.theta += param_1.dtheta;
	data1_1.x13 = *rtrack * cos(param_1.theta);
	data1_1.y13 = msign * *rtrack * sin(param_1.theta);
	transbwd_(x0, y0, z0);
	chkcmpt_();
	index_1.mroll = 2;
	findxyz_();
	if (index_1.movlpcond == 1) {
	    goto L10;
	}
	if (data_1.z1 < param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	} else if (data_1.x1 < param_1.bound || data_1.x1 > param_1.size - 
		param_1.bound || data_1.y1 < param_1.bound || data_1.y1 > 
		param_1.size - param_1.bound) {
	    index_1.nrec = 1;
	    goto L20;
	} else if (index_1.mroll == 3) {
	    if (param_1.dtheta <= param_1.eps) {
		data_1.x1 = x1tmp;
		data_1.y1 = y1tmp;
		data_1.z1 = z1tmp;
		data_1.ix1 = ix1tmp;
		data_1.jy1 = jy1tmp;
		data_1.kz1 = kz1tmp;
	    }
	    goto L20;
	}
    }
    index_1.mroll = 1;
    if (data_1.z3 < data_1.z2) {
	data_1.x2 = data_1.x3;
	data_1.y2 = data_1.y3;
	data_1.z2 = data_1.z3;
	data_1.r2 = data_1.r3;
	data_1.ix2 = data_1.ix3;
	data_1.jy2 = data_1.jy3;
	data_1.kz2 = data_1.kz3;
	data_1.icst2 = data_1.icst3;
    }
L20:
    return 0;
} /* calrloc_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int transfwd_(real *x1max, real *y1max, real *z1max, real *
	x0, real *y0, real *z0, real *rtrack)
{
    /* Builtin functions */
    double sqrt(doublereal);

/* .....Parallel translation */
    data1_1.x11max = *x1max - *x0;
    data1_1.y11max = *y1max - *y0;
    data1_1.z11max = *z1max - *z0;
    data1_1.x11 = data_1.x1 - *x0;
    data1_1.y11 = data_1.y1 - *y0;
    data1_1.z11 = data_1.z1 - *z0;
/* .....Rotation of x and y axes by angle A */
    data1_1.rxy = sqrt(*rtrack * *rtrack - data1_1.z11max * data1_1.z11max);
    data1_1.sinea = data1_1.y11max / data1_1.rxy;
    if (data1_1.sinea > 1.f) {
	data1_1.sinea = 1.f;
    }
    if (data1_1.sinea < -1.f) {
	data1_1.sinea = -1.f;
    }
    if (data1_1.x11max < 0.f) {
	data1_1.cosinea = -sqrt(1.f - data1_1.sinea * data1_1.sinea);
    } else {
	data1_1.cosinea = sqrt(1.f - data1_1.sinea * data1_1.sinea);
    }
    data1_1.x12 = data1_1.x11 * data1_1.cosinea + data1_1.y11 * data1_1.sinea;
    data1_1.y12 = -data1_1.x11 * data1_1.sinea + data1_1.y11 * 
	    data1_1.cosinea;
    data1_1.z12 = data1_1.z11;
/* .....Rotation of x and z axes by angle B */
    data1_1.sineb = data1_1.z11max / *rtrack;
    if (data1_1.sineb > 1.f) {
	data1_1.sineb = 1.f;
    }
    data1_1.cosineb = sqrt(1.f - data1_1.sineb * data1_1.sineb);
    data1_1.x13 = data1_1.x12 * data1_1.cosineb + data1_1.z12 * data1_1.sineb;
    data1_1.y13 = data1_1.y12;
    return 0;
} /* transfwd_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int transbwd_(real *x0, real *y0, real *z0)
{
/* .....Rotation of x and z axes by angle B */
    data1_1.x12 = data1_1.x13 * data1_1.cosineb;
    data1_1.y12 = data1_1.y13;
    data1_1.z12 = data1_1.x13 * data1_1.sineb;
/* .....Rotation of x and y axes by angle A */
    data1_1.x11 = data1_1.x12 * data1_1.cosinea - data1_1.y12 * data1_1.sinea;
    data1_1.y11 = data1_1.x12 * data1_1.sinea + data1_1.y12 * data1_1.cosinea;
    data1_1.z11 = data1_1.z12;
/* .....Parallel translation */
    data_1.x1 = data1_1.x11 + *x0;
    data_1.y1 = data1_1.y11 + *y0;
    data_1.z1 = data1_1.z11 + *z0;
    return 0;
} /* transbwd_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int chkcond_(void)
{
    /* Builtin functions */
    double sqrt(doublereal), acos(doublereal);

    /* Local variables */
    static real d12, d13, d14, d23, d24, d34, case1, case2, case3, alpha1, 
	    alpha2, alpha3, alpha4, ratio1, ratio2, ratio3, ratio4, alpha12, 
	    alpha34, ratio12, ratio34;

    d12 = sqrt((data_1.x1 - data_1.x2) * (data_1.x1 - data_1.x2) + (data_1.y1 
	    - data_1.y2) * (data_1.y1 - data_1.y2));
    d13 = sqrt((data_1.x1 - data_1.x3) * (data_1.x1 - data_1.x3) + (data_1.y1 
	    - data_1.y3) * (data_1.y1 - data_1.y3));
    d14 = sqrt((data_1.x1 - data_1.x4) * (data_1.x1 - data_1.x4) + (data_1.y1 
	    - data_1.y4) * (data_1.y1 - data_1.y4));
    d23 = sqrt((data_1.x2 - data_1.x3) * (data_1.x2 - data_1.x3) + (data_1.y2 
	    - data_1.y3) * (data_1.y2 - data_1.y3));
    d24 = sqrt((data_1.x2 - data_1.x4) * (data_1.x2 - data_1.x4) + (data_1.y2 
	    - data_1.y4) * (data_1.y2 - data_1.y4));
    d34 = sqrt((data_1.x3 - data_1.x4) * (data_1.x3 - data_1.x4) + (data_1.y3 
	    - data_1.y4) * (data_1.y3 - data_1.y4));
    ratio1 = (d23 * d23 + d12 * d12 - d13 * d13) * .5f / (d23 * d12);
    ratio2 = (d24 * d24 + d12 * d12 - d14 * d14) * .5f / (d24 * d12);
    ratio3 = (d23 * d23 + d13 * d13 - d12 * d12) * .5f / (d23 * d13);
    ratio4 = (d34 * d34 + d13 * d13 - d14 * d14) * .5f / (d34 * d13);
    ratio12 = (d23 * d23 + d24 * d24 - d34 * d34) * .5f / (d23 * d24);
    ratio34 = (d23 * d23 + d34 * d34 - d24 * d24) * .5f / (d23 * d34);
    if (ratio1 > 1.f) {
	ratio1 = 1.f;
    }
    if (ratio2 > 1.f) {
	ratio2 = 1.f;
    }
    if (ratio3 > 1.f) {
	ratio3 = 1.f;
    }
    if (ratio4 > 1.f) {
	ratio4 = 1.f;
    }
    if (ratio12 > 1.f) {
	ratio12 = 1.f;
    }
    if (ratio34 > 1.f) {
	ratio34 = 1.f;
    }
    if (d14 <= param_1.tolrnce) {
	index_1.nrec = 1;
    } else {
	alpha1 = acos(ratio1);
	alpha2 = acos(ratio2);
	alpha3 = acos(ratio3);
	alpha4 = acos(ratio4);
	alpha12 = acos(ratio12);
	alpha34 = acos(ratio34);
	case1 = alpha2 - alpha12;
	case2 = alpha1 - alpha12;
	case3 = alpha3 - alpha34;
	if (case1 > param_1.tolrnce) {
	    data_1.x1 += 1e-4f;
	    if (data_1.x1 >= param_1.cmptln * data_1.ix1 && data_1.ix1 < 
		    param_1.nc) {
		++data_1.ix1;
	    }
	    data_1.y1 += 1e-4f;
	    if (data_1.y1 >= param_1.cmptln * data_1.jy1 && data_1.jy1 < 
		    param_1.nc) {
		++data_1.jy1;
	    }
	    index_1.mroll = 0;
	} else if (case2 > param_1.tolrnce) {
	    data_1.x3 = data_1.x4;
	    data_1.y3 = data_1.y4;
	    data_1.z3 = data_1.z4;
	    data_1.r3 = data_1.r4;
	    data_1.ix3 = data_1.ix4;
	    data_1.jy3 = data_1.jy4;
	    data_1.kz3 = data_1.kz4;
	    data_1.icst3 = data_1.icst4;
	    index_1.mroll = 2;
	} else if (case3 > param_1.tolrnce) {
	    data_1.x2 = data_1.x4;
	    data_1.y2 = data_1.y4;
	    data_1.z2 = data_1.z4;
	    data_1.r2 = data_1.r4;
	    data_1.ix2 = data_1.ix4;
	    data_1.jy2 = data_1.jy4;
	    data_1.kz2 = data_1.kz4;
	    data_1.icst2 = data_1.icst4;
	    index_1.mroll = 2;
	} else {
	    index_1.nrec = 1;
	}
    }
    return 0;
} /* chkcond_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int tmprec_(void)
{
    extern /* Subroutine */ int localshift_(void), choosesite_(void);

    local_1.move = 10;
    if (local_1.intensity == 0) {
	goto L10;
    }
    if (local_1.lshift <= local_1.move) {
	localshift_();
    }
    if (local_1.lshift == local_1.move) {
	choosesite_();
	goto L10;
    }
    ++local_1.lshift;
L10:
    return 0;
} /* tmprec_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int localshift_(void)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    extern doublereal ran_(integer *);
    static real sample;
    extern /* Subroutine */ int chkcmpt_(void);

    shift_1.xshift[local_1.lshift - 1] = data_1.x1;
    shift_1.yshift[local_1.lshift - 1] = data_1.y1;
    shift_1.zshift[local_1.lshift - 1] = data_1.z1;
    shift_1.ixshift[local_1.lshift - 1] = data_1.ix1;
    shift_1.jyshift[local_1.lshift - 1] = data_1.jy1;
    shift_1.kzshift[local_1.lshift - 1] = data_1.kz1;
    if (local_1.lshift < local_1.move) {
	sample = param_1.size;
	while(sample >= param_1.size - param_1.bound || sample <= 
		param_1.bound) {
	    sample = data_1.x1 - dist_1.rmax + (r__1 = dist_1.rmax * 2.f * 
		    ran_(&param_1.nrndsd), dabs(r__1));
	}
	data_1.x1 = sample;
	sample = param_1.size;
	while(sample >= param_1.size - param_1.bound || sample <= 
		param_1.bound) {
	    sample = data_1.y1 - dist_1.rmax + (r__1 = dist_1.rmax * 2.f * 
		    ran_(&param_1.nrndsd), dabs(r__1));
	}
	data_1.y1 = sample;
	chkcmpt_();
	index_1.mroll = 0;
	index_1.nrec = 0;
    }
    return 0;
} /* localshift_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int choosesite_(void)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, l;
    static real q;
    static integer ir, indxt;

    i__1 = local_1.move;
    for (j = 1; j <= i__1; ++j) {
	shift_1.indx[j - 1] = j;
    }
    l = local_1.move / 2 + 1;
    ir = local_1.move;
L10:
    if (l > 1) {
	--l;
	indxt = shift_1.indx[l - 1];
	q = shift_1.zshift[indxt - 1];
    } else {
	indxt = shift_1.indx[ir - 1];
	q = shift_1.zshift[indxt - 1];
	shift_1.indx[ir - 1] = shift_1.indx[0];
	--ir;
	if (ir == 1) {
	    shift_1.indx[0] = indxt;
	    goto L30;
	}
    }
    i__ = l;
    j = l + l;
L20:
    if (j <= ir) {
	if (j < ir) {
	    if (shift_1.zshift[shift_1.indx[j - 1] - 1] < shift_1.zshift[
		    shift_1.indx[j] - 1]) {
		++j;
	    }
	}
	if (q < shift_1.zshift[shift_1.indx[j - 1] - 1]) {
	    shift_1.indx[i__ - 1] = shift_1.indx[j - 1];
	    i__ = j;
	    j += j;
	} else {
	    j = ir + 1;
	}
	goto L20;
    }
    shift_1.indx[i__ - 1] = indxt;
    goto L10;
L30:
    i__ = local_1.move + 1 - local_1.intensity;
    data_1.x1 = shift_1.xshift[shift_1.indx[i__ - 1] - 1];
    data_1.y1 = shift_1.yshift[shift_1.indx[i__ - 1] - 1];
    data_1.z1 = shift_1.zshift[shift_1.indx[i__ - 1] - 1];
    data_1.ix1 = shift_1.ixshift[shift_1.indx[i__ - 1] - 1];
    data_1.jy1 = shift_1.jyshift[shift_1.indx[i__ - 1] - 1];
    data_1.kz1 = shift_1.kzshift[shift_1.indx[i__ - 1] - 1];
    index_1.nrec = 1;
    return 0;
} /* choosesite_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int relocate_(void)
{
    static integer nchalf;

    nchalf = (integer) (param_1.nc * .5f);
    if (index_1.mcount == 40 || index_1.ntshift >= 60) {
	data_1.r1 *= .9f;
	index_1.mcount = 0;
	index_1.ntshift = 0;
    }
L10:
    if (index_1.mcount == 0) {
	if (data_1.ix1 <= nchalf) {
	    sign_1.lsign = 1;
	} else {
	    sign_1.lsign = -1;
	}
    }
    if (data_1.jy1 < param_1.nc) {
	data_1.y1 += param_1.cmptln;
	++data_1.jy1;
    } else {
	data_1.x1 += param_1.cmptln * sign_1.lsign;
	data_1.ix1 += sign_1.lsign;
	data_1.y1 = data_1.y1 + param_1.cmptln - param_1.size;
	data_1.jy1 = 1;
    }
    if (data_1.ix1 == 0 || data_1.ix1 == param_1.nc + 1) {
	data_1.x1 -= param_1.cmptln * 2.f * (real) sign_1.lsign;
	data_1.ix1 -= sign_1.lsign << 1;
	data_1.r1 *= .9f;
	index_1.mcount = 0;
	index_1.ntshift = 0;
	goto L10;
    }
    ++index_1.mcount;
    index_1.mroll = 0;
    index_1.nrec = 0;
    ++index_1.ntshift;
    return 0;
} /* relocate_ */

/* ------------------------------------------------------------------------ */
/* Subroutine */ int finalset_(void)
{
    /* Format strings */
    static char fmt_2000[] = "(\002                     Porosity of the samp"
	    "le: \002,e13.5,/,\002                              Sample length"
	    ": \002,e13.5,/,\002                               Sample width:"
	    " \002,e13.5,/,\002                               Sample depth:"
	    " \002,e13.5,/,\002                                Cell length:"
	    " \002,e13.5,/,\002                        Number of particles:"
	    " \002,i13,/,\002       Number of cells in each x and y axis: "
	    "\002,i13,/,\002                  Number of cells in z axis: \002"
	    ",i13,/,\002                Periodical BC (Yes=1, No=0): \002,i13)"
	    ;
    static char fmt_2010[] = "(\002                           Average of ln("
	    "r): \002,e13.5,/,\002                             S.D. for ln(r)"
	    ": \002,e13.5,/,\002                        Maximum samplable r:"
	    " \002,e13.5,/,\002                        Minimum samplable r:"
	    " \002,e13.5,/,\002                        Average samplable r:"
	    " \002,e13.5,/,\002                Small particle radius limit:"
	    " \002,e13.5)";
    static char fmt_2015[] = "(\002                        Number of iterati"
	    "on: \002,i13,/,\002                             Porosity index:"
	    " \002,i13)";
    static char fmt_1000[] = "(/,4i7)";
    static char fmt_2020[] = "(4i7)";
    static char fmt_2030[] = "(4(e13.5,2x))";
    static char fmt_2016[] = "(\002              Number of particles extract"
	    "ed: \002,i13,/,\002        Starting cell (in x or y direction):"
	    " \002,i13,/,\002          Ending cell (in x or y direction): "
	    "\002,i13,/,\002             Starting cell (in z direction): \002"
	    ",i13,/,\002               Ending cell (in z direction): \002,i13)"
	    ;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    olist o__1;
    alist al__1;

    /* Builtin functions */
    double exp(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     f_open(olist *), s_rsfe(cilist *), e_rsfe(void), f_rew(alist *);

    /* Local variables */
    static real r1, x1, y1, z1;
    static integer ic, ix, jy, kz, iic, nch, ncl;
    static real ravg;
    static integer icst, nczh, nczl;
    static real vpack, avgvc, hctop, trueh, vctop;
    static integer ncntst;

    /* Fortran I/O blocks */
    static cilist io___299 = { 0, 7, 0, fmt_2000, 0 };
    static cilist io___300 = { 0, 7, 0, fmt_2010, 0 };
    static cilist io___301 = { 0, 7, 0, fmt_2015, 0 };
    static cilist io___302 = { 0, 5, 0, fmt_1000, 0 };
    static cilist io___311 = { 0, 9, 0, fmt_2020, 0 };
    static cilist io___313 = { 0, 9, 0, fmt_2030, 0 };
    static cilist io___314 = { 0, 7, 0, fmt_2016, 0 };


/* .....Calculate the final porosity, number of spheres in the whole */
/*     domain */
    vpack = 0.f;
    param_1.npart = 0;
    ravg = exp(dist_1.alphar);
    dist_1.rmin = ravg;
    dist_1.rmax = ravg;
    i__1 = param_1.ncz - 2;
    for (kz = 1; kz <= i__1; ++kz) {
	i__2 = param_1.nc;
	for (jy = 1; jy <= i__2; ++jy) {
	    i__3 = param_1.nc;
	    for (ix = 1; ix <= i__3; ++ix) {
		ncntst = sortdim_1.icntst[ix + (jy + kz * 40) * 40 - 1641];
		i__4 = ncntst;
		for (icst = 1; icst <= i__4; ++icst) {
		    r1 = sortdim_1.rcst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    vpack += param_1.pi * 4.f * r1 * r1 * r1 / 3.f;
		    if (r1 < dist_1.rmin) {
			dist_1.rmin = r1;
		    }
		    if (r1 > dist_1.rmax) {
			dist_1.rmax = r1;
		    }
		    ++param_1.npart;
		}
	    }
	}
    }
    avgvc = vpack / (param_1.ncz - 2.f);
    vctop = 0.f;
    i__1 = param_1.ncz;
    for (kz = param_1.ncz - 1; kz <= i__1; ++kz) {
	i__2 = param_1.nc;
	for (jy = 1; jy <= i__2; ++jy) {
	    i__3 = param_1.nc;
	    for (ix = 1; ix <= i__3; ++ix) {
		ncntst = sortdim_1.icntst[ix + (jy + kz * 40) * 40 - 1641];
		i__4 = ncntst;
		for (icst = 1; icst <= i__4; ++icst) {
		    r1 = sortdim_1.rcst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    vctop += param_1.pi * 4.f * r1 * r1 * r1 / 3.f;
		    if (r1 < dist_1.rmin) {
			dist_1.rmin = r1;
		    }
		    if (r1 > dist_1.rmax) {
			dist_1.rmax = r1;
		    }
		    ++param_1.npart;
		}
	    }
	}
    }
    hctop = param_1.cmptln * vctop / avgvc;
    if (hctop >= param_1.cmptln * 2) {
	hctop = param_1.cmptln * 2;
    }
    trueh = (param_1.ncz - 2) * param_1.cmptln + hctop;
    param_1.porosity = 1.f - avgvc / (param_1.size * param_1.size * 
	    param_1.cmptln);
    s_wsfe(&io___299);
    do_fio(&c__1, (char *)&param_1.porosity, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.size, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.size, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&trueh, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.cmptln, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.npart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&param_1.nc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&param_1.ncz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&index_1.mpbc, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___300);
    do_fio(&c__1, (char *)&dist_1.alphar, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&dist_1.betar, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&dist_1.rmax, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&dist_1.rmin, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&ravg, (ftnlen)sizeof(real));
    do_fio(&c__1, (char *)&param_1.rlmt, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___301);
    do_fio(&c__1, (char *)&index_1.niterate, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&local_1.intensity, (ftnlen)sizeof(integer));
    e_wsfe();
/* .....Write out the final extracted packing data! */
    o__1.oerr = 0;
    o__1.ounit = 5;
    o__1.ofnmlen = 10;
    o__1.ofnm = "subdim.par";
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsfe(&io___302);
    do_fio(&c__1, (char *)&ncl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nczl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nczh, (ftnlen)sizeof(integer));
    e_rsfe();
    al__1.aerr = 0;
    al__1.aunit = 9;
    f_rew(&al__1);
    param_1.npart = 0;
    i__1 = nczh;
    for (kz = nczl; kz <= i__1; ++kz) {
	i__2 = nch;
	for (jy = ncl; jy <= i__2; ++jy) {
	    i__3 = nch;
	    for (ix = ncl; ix <= i__3; ++ix) {
		ncntst = sortdim_1.icntst[ix + (jy + kz * 40) * 40 - 1641];
		ic = 0;
		i__4 = ncntst;
		for (icst = 1; icst <= i__4; ++icst) {
		    x1 = sortdim_1.xcst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    y1 = sortdim_1.ycst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    z1 = sortdim_1.zcst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    r1 = sortdim_1.rcst[ix + (jy + (kz + icst * 40) * 40) * 
			    40 - 65641];
		    if (x1 >= param_1.tolrnce && x1 <= param_1.size + 
			    param_1.tolrnce && y1 >= param_1.tolrnce && y1 <= 
			    param_1.size + param_1.tolrnce && z1 >= 
			    param_1.tolrnce && z1 <= param_1.height + 
			    param_1.tolrnce) {
			++param_1.npart;
			++ic;
			origdim_1.x[ic - 1] = x1;
			origdim_1.y[ic - 1] = y1;
			origdim_1.z__[ic - 1] = z1;
			origdim_1.r__[ic - 1] = r1;
		    }
		}
		s_wsfe(&io___311);
		do_fio(&c__1, (char *)&ix, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&jy, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&kz, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ic, (ftnlen)sizeof(integer));
		e_wsfe();
		i__4 = ic;
		for (iic = 1; iic <= i__4; ++iic) {
		    s_wsfe(&io___313);
		    do_fio(&c__1, (char *)&origdim_1.x[iic - 1], (ftnlen)
			    sizeof(real));
		    do_fio(&c__1, (char *)&origdim_1.y[iic - 1], (ftnlen)
			    sizeof(real));
		    do_fio(&c__1, (char *)&origdim_1.z__[iic - 1], (ftnlen)
			    sizeof(real));
		    do_fio(&c__1, (char *)&origdim_1.r__[iic - 1], (ftnlen)
			    sizeof(real));
		    e_wsfe();
		}
	    }
	}
    }
    s_wsfe(&io___314);
    do_fio(&c__1, (char *)&param_1.npart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ncl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nch, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nczl, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nczh, (ftnlen)sizeof(integer));
    e_wsfe();
    return 0;
} /* finalset_ */

/* Main program alias */ int randompack_ () { MAIN__ (); return 0; }
