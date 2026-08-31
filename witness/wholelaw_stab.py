"""wholelaw.py + DUAL STABILISATION.  Drop-in replacement; same CLI:

    python wholelaw_stab.py <tag> <outers> <rounds> <reps> <nproc> [alpha]

Two changes over wholelaw.py, both measured (see the dual-stabilisation report):

 1. PRICING DUAL = Wentges / in-out smoothing.  The pricing dual is
        q_price = alpha*q_LP + (1-alpha)*qbar,     qbar <- 0.5*qbar + 0.5*q_LP
    with alpha = 0.5 by default.  Measured on an identical 3-round A/B from the wl3
    state at obj 0.681955776635: total decrease 3.719e-06 vs 3.093e-06 for the plain
    LP dual (+20%), with the extra gain coming from the LP being able to put MORE weight
    on the newly generated columns (w_new 0.365/0.262 vs 0.264/0.218), not from deeper
    reduced costs.
    IMPORTANT: the smoothed dual must stay ON the manifold of exact duals of nearby
    masters.  An ISOTROPIC RANDOM perturbation of the same size (|dq|_inf = 2e-7) kills
    pricing outright: 0 of 256 slid columns price negative and the objective does not
    move at all.  Do not substitute "add noise to q" for this.

 2. The LP dual is cross-checked against the exact dual obtained from complementary
    slackness (a 256x256 solve).  The pricing oracle tolerates |dq|_inf <~ 1e-8; HiGHS
    at its DEFAULT tolerances returns 8.1e-8, which is already destructive.  The
    dual_feasibility_tolerance=1e-10 option in the master is load-bearing.  Note that
    asking for 1e-12 is rejected by HiGHS as an invalid value and SILENTLY falls back to
    the defaults -- do not "tighten" it.
"""
import os,sys,time,pickle
os.environ.setdefault("OMP_NUM_THREADS","1"); os.environ.setdefault("MKL_NUM_THREADS","1")
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog
from colgen import ROWS, rows_of_config, sigma_of
from sfw import slide
from duals import extract
from exactdual import extract_exact
R=ROWS.astype(float)
DUAL_BUDGET=1e-8
Q=None;MU=None;LAM=0.0
def _init(q,mu,lam=0.0):
    # RE-APPLIED AUDIT FIX (bug 15): lam, the dual of the F(256) floor row, must reach the
    # workers.  The master prices sigma + q.F + mu + lam*F(256); without this the slide
    # optimises a different objective whenever that row binds.  Harmless only while lam == 0.
    global Q,MU,LAM; Q=q;MU=mu;LAM=lam
def _job(a):
    i,x,m,rounds,jit,seed=a
    rng=np.random.default_rng(seed)
    x=np.asarray(x,float).copy(); m=np.asarray(m,float)
    if jit>0: x=np.mod(x+rng.normal(0,jit,len(x)),256.0)
    best=None
    for r in range(rounds):
        x,E=slide(x,m,Q,MU); rc=sigma_of(m)+E+MU
        if LAM: rc += LAM*f2(x,m)
        if best is None or rc<best[0]: best=(rc,x.copy())
    return (i,float(best[0]),best[1],m)
def f2(x,m):
    a=2*np.pi*np.asarray(x,float); return float((np.cos(a)*m).sum()**2+(np.sin(a)*m).sum()**2)

if __name__=="__main__":
    tag=sys.argv[1]; outers=int(sys.argv[2]); rounds=int(sys.argv[3]); reps=int(sys.argv[4]); nproc=int(sys.argv[5])
    ALPHA=float(sys.argv[6]) if len(sys.argv)>6 else 0.5
    A=np.load(f"{tag}_A.npy"); cost=np.load(f"{tag}_cost.npy"); fv=np.load(f"{tag}_f256.npy")
    POS=pickle.load(open(f"{tag}_pos.pkl","rb"))
    sup=np.load(f"{tag}_supidx.npy"); q=np.load(f"{tag}_q.npy")
    dd=pickle.load(open(f"{tag}_dual.pkl","rb")); mu=dd["mu"]; obj=dd["obj"]
    # D(1) is a CONSTRAINT of the theorem, not a post-hoc check.  E[F(256)] is a free direction
    # for lowering p, so an uncapped solve WILL drift above the published cap and produce a law
    # that is inadmissible -- discovered only after an hour of interval arithmetic.  Carry both
    # sides here.  Costs about 4.6e-7 in objective; measured.
    from fractions import Fraction as _Fr
    F256_HI=float(_Fr(128)+_Fr(82395317,100000000)*65536)   # E[F(256)] <= 54126.59494912
    b=np.concatenate([R+1e-4,-(R-1e-4),[-128.0],[F256_HI]])
    active=np.union1d(sup, np.argsort(cost+q@A+mu)[:2500])
    qbar=q.copy(); mubar=mu
    t0=time.time()
    def solve(ix):
        Au=A[:,ix]; f=fv[ix]
        return linprog(cost[ix],A_ub=np.vstack([Au,-Au,-f[None,:],f[None,:]]),b_ub=b,A_eq=np.ones((1,len(ix))),b_eq=[1.0],
            bounds=(0,None),method="highs",options={"dual_feasibility_tolerance":1e-10,"primal_feasibility_tolerance":1e-10}),Au
    lam0=float(dd.get("lam",0.0))
    for outer in range(outers):
        qp = ALPHA*q + (1.0-ALPHA)*qbar          # <-- stabilised PRICING dual
        mup = ALPHA*mu + (1.0-ALPHA)*mubar
        jobs=[]
        for k,ci in enumerate(sup):
            if POS.get(int(ci)) is None: continue
            x,m=POS[int(ci)]
            for rp in range(reps):
                jobs.append((k,x,m,rounds,0.0 if rp==0 else 0.015*rp,hash((outer,k,rp))%(1<<62)))
        res={}
        with Pool(nproc,initializer=_init,initargs=(qp,mup,lam0)) as P:
            for i,rc,x,m in P.imap_unordered(_job,jobs,chunksize=2):
                if i not in res or rc<res[i][0]: res[i]=(rc,x,m)
        ks=sorted(res)
        rcs=np.array([res[k][0] for k in ks])
        newA=np.array([rows_of_config(res[k][1],res[k][2]) for k in ks]).T
        newc=np.array([sigma_of(res[k][2]) for k in ks])
        newf=np.array([f2(res[k][1],res[k][2]) for k in ks])
        base=A.shape[1]
        A=np.hstack([A,newA]); cost=np.concatenate([cost,newc]); fv=np.concatenate([fv,newf])
        for n_,k in enumerate(ks): POS[base+n_]=(res[k][1],res[k][2])
        active=np.union1d(active,np.arange(base,base+len(ks)))
        rc_true=newc+q@newA+mu+lam0*newf
        print(f"outer {outer}: slid {len(ks)}  |dq_smooth| {np.abs(qp-q).max():.2e}  "
              f"rc(exact dual) min {rc_true.min():+.3e} med {np.median(rc_true):+.3e} "
              f"neg {int((rc_true<0).sum())}/{len(ks)}  [{time.time()-t0:.0f}s]",flush=True)
        for it in range(30):
            r,Au=solve(active)
            if not r.success: print("  INFEASIBLE"); sys.exit(1)
            q,mu,rep=extract(r,Au,cost[active],f256=fv[active])
            rc=cost+q@A+mu+rep["lam"]*fv
            neg=np.where(rc<-1e-12)[0]
            if len(neg)==0: break
            add=np.setdiff1d(neg[np.argsort(rc[neg])][:2500],active)
            if len(add)==0: break
            supn=active[np.where(r.x>1e-11)[0]]
            active=np.union1d(np.union1d(supn,add),active[np.argsort(rc[active])[:2500]])
        obj=float(r.fun); supl=np.where(r.x>1e-11)[0]; sup=active[supl]
        # --- dual accuracy audit: pricing needs |dq|_inf <~ 1e-8 ---
        try:
            qe,mue,info=extract_exact(A,cost,sup,f256=fv,lam=rep["lam"])
            d=float(np.abs(qe-q).max())
            if d>DUAL_BUDGET:
                print(f"   !! DUAL ERROR {d:.2e} EXCEEDS THE {DUAL_BUDGET:.0e} PRICING BUDGET "
                      f"(CS residual {info['residual']:.1e}); using the exact dual",flush=True)
                q,mu=qe,mue
            else:
                print(f"   dual audit ok: |q_LP - q_exact| = {d:.2e}  (budget {DUAL_BUDGET:.0e}, "
                      f"CS residual {info['residual']:.1e}, cond {info['cond']:.1e})",flush=True)
        except RuntimeError as e:
            print(f"   dual audit skipped: {e}",flush=True)
        qbar=0.5*qbar+0.5*q; mubar=0.5*mubar+0.5*mu
        print(f"   -> obj {obj:.12f}   gap to p0 {obj-0.6818286874638315:+.3e}   support {len(sup)}  "
              f"w_new {float(r.x[supl][sup>=base].sum()):.4f}  [{time.time()-t0:.0f}s]",flush=True)
        np.save(f"{tag}_A.npy",A); np.save(f"{tag}_cost.npy",cost); np.save(f"{tag}_f256.npy",fv)
        np.save(f"{tag}_supidx.npy",sup); np.save(f"{tag}_w.npy",r.x[supl]); np.save(f"{tag}_q.npy",q)
        lam0=float(rep["lam"])
        pickle.dump({"mu":mu,"lam":rep["lam"],"obj":obj},open(f"{tag}_dual.pkl","wb"))
        pickle.dump(POS,open(f"{tag}_pos.pkl","wb"),protocol=4)
