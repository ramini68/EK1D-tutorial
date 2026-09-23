// Optional standard-library C++17 implementation of the corrected Python solver.
// No fast-math, concentration clipping, or spectator-ion charge adjustment.
// c is species-major H, OH, Cu, Na, NO3; all states are mol/m3 pore water.
// params: L,theta,tau,V,Kw,Ksp,pH50,nedge,dt,cfl,D0[5],initial_inventory[3],spatial_order,algebraic_n2,time_order.
// metrics: steps,min_dt,max_dt,min_c,max_abs_charge,max_charge_fraction,
// max_current_residual,max_voltage_residual,relative_inventory_error[3],
// water_residual,sorption_residual,precipitation_residual,Faradaic_moles,
// advanced_duration,max_local_charge_residual,max_local_Cu_residual,max_chem_iter,
// algebraic_chemistry_successes,generic_chemistry_calls,rejected_RK2_attempts,
// accepted_stage_current_min,accepted_stage_current_max,signed_charge_C_m2.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <utility>

namespace {
constexpr double F=96485.0, FRT=96485.0/(8.314*298.15);
constexpr double z[5]={1.0,-1.0,2.0,1.0,-1.0};
constexpr double eps=std::numeric_limits<double>::epsilon();
double water_h(double d,double kw) {
    const double r=std::hypot(d,2.0*std::sqrt(kw));
    return d>=0.0 ? .5*d+.5*r : kw/(.5*r-.5*d);
}
double minmod_mc(double a,double b,double c) {
    if(a>0.0&&b>0.0&&c>0.0)return std::min({a,b,c});
    if(a<0.0&&b<0.0&&c<0.0)return std::max({a,b,c});
    return 0.0;
}
struct Chem {
    double kw,ksp,n,lh50,lkw,lksp,h50,kwc,h50c,B,cs_sat;
    bool algebraic;
    explicit Chem(const double* p):kw(p[4]),ksp(p[5]),n(p[7]),
        lh50(-p[6]*std::log(10.0)),lkw(std::log(p[4])),lksp(std::log(p[5])),
        h50(std::pow(10.0,-p[6])),kwc(p[4]*1e6),h50c(h50*1000.0),
        B((p[5]*1e9)/(kwc*kwc)),cs_sat(B*h50c*h50c),algebraic(p[19]!=0.0) {}
    bool solve_n2(double total,double charge,double hint,double& h_out,double& oh_out,
                  double& aq_out,double& sorb_out,double& ppt_out,double* metrics) const {
        // Algebraic evaluation of precisely the same n=2 scalar equilibrium.
        // Calculations use mol/m3 throughout; no reaction is approximated.
        if(!(kwc>0.0&&h50c>0.0&&B>0.0)||!std::isfinite(cs_sat))return false;
        double lo=water_h(charge-2.0*total,kwc),hi=water_h(charge,kwc);
        double h=std::sqrt(lo)*std::sqrt(hi),previous=hi-lo;
        if(std::isfinite(hint)&&hint>0.0)h=std::min(hi,std::max(lo,hint));
        for(int it=1;it<=40;++it) {
            if(!(h>0.0)||!std::isfinite(h))return false;
            const double oh=kwc/h;
            const double ratio=h>=h50c?h50c/h:h/h50c,r2=ratio*ratio;
            const double af=h>=h50c?1.0/(1.0+r2):r2/(1.0+r2);
            const double sf=h>=h50c?r2/(1.0+r2):1.0/(1.0+r2);
            const double au=total*af,su=total*sf,as=B*h*h;
            bool saturated=total>0.0&&as<au;
            double aq=saturated?as:au,sorb=saturated?cs_sat:su;
            double ppt=saturated?total-(aq+sorb):0.0;
            if(ppt<0.0){if(ppt < -64.0*eps*total)return false;aq=au;sorb=su;ppt=0.0;saturated=false;}
            const double res=(h-oh)+2.0*aq-charge;
            const double scale=h+oh+2.0*aq+std::abs(charge);
            if(!std::isfinite(res)||!std::isfinite(scale))return false;
            if(std::abs(res)<=64.0*eps*scale+1e-21) {
                const double cr=aq+sorb+ppt-total;
                if(std::abs(cr)>64.0*eps*total+1e-15)return false;
                h_out=h;oh_out=oh;aq_out=aq;sorb_out=sorb;ppt_out=ppt;
                metrics[16]=std::max(metrics[16],std::abs(res));metrics[17]=std::max(metrics[17],std::abs(cr));
                metrics[18]=std::max(metrics[18],double(it));metrics[19]+=1.0;
                return true;
            }
            if(res<0.0)lo=h;
            if(res>0.0)hi=h;
            const double deriv=1.0+(oh+4.0*aq*(saturated?1.0:sf))/h;
            double proposal=h-res/deriv;
            if(!std::isfinite(proposal)||proposal<=lo||proposal>=hi||std::abs(proposal-h)>.5*previous)
                proposal=.5*lo+.5*hi;
            previous=std::abs(proposal-h);h=proposal;
        }
        return false;
    }
    void solve(double total,double charge,double hint,double& h_out,double& oh_out,
               double& aq_out,double& sorb_out,double& ppt_out,double* metrics) const {
        if(algebraic&&n==2.0&&solve_n2(total,charge,hint,h_out,oh_out,aq_out,sorb_out,ppt_out,metrics))return;
        metrics[20]+=1.0;
        const double t=total/1000.0,a=charge/1000.0;
        const double logt=t>0.0 ? std::log(t) : -std::numeric_limits<double>::infinity();
        double lo=std::log(water_h(a-2.0*t,kw)),hi=std::log(water_h(a,kw));
        double y=.5*(lo+hi),previous=hi-lo;
        if(std::isfinite(hint)&&hint>0.0) y=std::min(hi,std::max(lo,std::log(hint)-std::log(1000.0)));
        double h=0,oh=0,aq=0,sorb=0,ppt=0,res=0;
        int it;
        for(it=1;it<=80;++it) {
            h=std::exp(y); oh=std::exp(lkw-y);
            const double u=n*(y-lh50),ex=std::exp(-std::abs(u)),den=1.0+ex;
            const double af=u>=0.0 ? 1.0/den : ex/den;
            const double sf=u>=0.0 ? ex/den : 1.0/den;
            const double au=t*af,su=t*sf;
            const double las=lksp+2.0*y-2.0*lkw;
            const double lau=logt-(std::max(0.0,-u)+std::log1p(std::exp(-std::abs(u))));
            bool saturated=t>0.0&&las<lau;
            aq=au; sorb=su; ppt=0.0;
            if(saturated) {
                aq=std::exp(las); sorb=std::exp(las-u); ppt=t-(aq+sorb);
                if(ppt<0.0) {
                    if(ppt < -64.0*eps*t) throw std::runtime_error("Inconsistent phase selection");
                    aq=au;sorb=su;ppt=0.0;saturated=false;
                }
            }
            res=(h-oh)+2.0*aq-a;
            const double deriv=h+oh+(saturated ? 4.0*aq : 2.0*n*aq*sf);
            const double scale=h+oh+2.0*aq+std::abs(a);
            if(std::abs(res)<=64.0*eps*scale+1e-24) break;
            if(res<0.0) lo=y;
            if(res>0.0) hi=y;
            double proposal=y-res/deriv;
            if(!std::isfinite(proposal)||proposal<=lo||proposal>=hi||std::abs(proposal-y)>.5*previous)
                proposal=.5*(lo+hi);
            previous=std::abs(proposal-y);y=proposal;
        }
        if(it>80) {
            char msg[256];std::snprintf(msg,sizeof(msg),"Chemistry did not converge: T=%.17g A=%.17g H=%.17g residual=%.17g",total,charge,h*1000.0,res*1000.0);
            throw std::runtime_error(msg);
        }
        const double cr=aq+sorb+ppt-t;
        if(std::abs(cr)>64.0*eps*t+1e-18) throw std::runtime_error("Chemistry copper conservation failed");
        h_out=h*1000.0;oh_out=oh*1000.0;aq_out=aq*1000.0;sorb_out=sorb*1000.0;ppt_out=ppt*1000.0;
        metrics[16]=std::max(metrics[16],std::abs(res)*1000.0);
        metrics[17]=std::max(metrics[17],std::abs(cr)*1000.0);
        metrics[18]=std::max(metrics[18],double(it));
    }
};
}

extern "C" int advance(int N,double* c,double* cs,double* cp,double duration,
                       const double* p,double* m,char* err,int errlen) {
    double elapsed=0.0;
    try {
        if(N<2||!std::isfinite(duration)||duration<0.0||p[0]<=0.0||p[1]<=0.0||p[2]<=0.0||
           p[4]<=0.0||p[5]<=0.0||p[7]<=0.0||p[8]<=0.0||p[9]<=0.0||p[9]>=1.0)
            throw std::runtime_error("Invalid solver parameters");
        for(int k=0;k<21;++k)if(!std::isfinite(p[k]))throw std::runtime_error("Nonfinite solver parameter");
        if(p[18]!=1.0&&p[18]!=2.0)throw std::runtime_error("Spatial order must be 1 or 2");
        if(p[20]!=1.0&&p[20]!=2.0)throw std::runtime_error("Time order must be 1 or 2");
        for(int k=15;k<18;++k)if(p[k]<=0.0)throw std::runtime_error("Initial inventories must be positive");
        for(int k=0;k<5*N;++k)if(c[k]<0.0||!std::isfinite(c[k]))throw std::runtime_error("Invalid initial aqueous state");
        for(int k=0;k<N;++k)if(cs[k]<0.0||cp[k]<0.0||!std::isfinite(cs[k]+cp[k]))throw std::runtime_error("Invalid initial solid state");
        std::fill(m,m+25,0.0);m[1]=std::numeric_limits<double>::infinity();m[3]=*std::min_element(c,c+5*N);
        m[22]=std::numeric_limits<double>::infinity();m[23]=-std::numeric_limits<double>::infinity();
        const int nf=N-1;const double dx=p[0]/N,theta=p[1],scale=theta*dx;
        double D[5],factor[5];for(int s=0;s<5;++s){D[s]=p[10+s]*theta/p[2];factor[s]=F*FRT*z[s]*z[s]*D[s];}
        Chem chem(p);
        std::vector<double> sigma(N),sp(nf),sm(nf),diff(nf),left(5*nf),right(5*nf),flux(5*(N+1)),next(5*N),sf(nf),slope(5*N),recl(5*nf),recr(5*nf);
        std::vector<int> order(nf);
        const auto transport_proposal = [&](const double* c,double limit) {
            if(p[18]==2.0) {
                std::fill(slope.begin(),slope.end(),0.0);
                for(int s=0;s<5;++s)for(int i=1;i<N-1;++i){const int k=s*N+i;const double dl=c[k]-c[k-1],dr=c[k+1]-c[k];
                    slope[k]=minmod_mc(2.0*dl,.5*(dl+dr),2.0*dr);}
            }
            for(int s=0;s<5;++s)for(int f=0;f<nf;++f){const int k=s*nf+f;
                recl[k]=c[s*N+f]+.5*slope[s*N+f];recr[k]=c[s*N+f+1]-.5*slope[s*N+f+1];
                if(recl[k]<0.0||recr[k]<0.0)throw std::runtime_error("Negative reconstructed migration concentration");}
            for(int i=0;i<N;++i){double v=0;for(int s=0;s<5;++s)v+=factor[s]*c[s*N+i];sigma[i]=v;if(!(v>0.0))throw std::runtime_error("Nonpositive physical conductivity");}
            for(int f=0;f<nf;++f) {
                double a=0,b=0,d=0;
                for(int s=0;s<5;++s){const double l=c[s*N+f],r=c[s*N+f+1],rl=recl[s*nf+f],rrr=recr[s*nf+f];
                    a+=factor[s]*(z[s]>0?rl:rrr);b+=factor[s]*(z[s]>0?rrr:rl);d+=z[s]*D[s]*(r-l)/dx;}
                sp[f]=a;sm[f]=b;diff[f]=-F*d;
                if(!(a>0.0&&b>0.0))throw std::runtime_error("Nonpositive face conductivity");
            }
            std::iota(order.begin(),order.end(),0);
            std::sort(order.begin(),order.end(),[&](int a,int b){return diff[a]<diff[b];});
            const double rb=dx*.5*(1.0/sigma[0]+1.0/sigma[N-1]);
            double rr=rb,rd=0.0;for(int f:order){const double r=dx/sm[f];rr+=r;rd+=r*diff[f];}
            double current=0.0;bool found=false;
            for(int k=0;k<=nf;++k) {
                const double candidate=(p[3]+rd)/rr;
                const double lo=k==0 ? -std::numeric_limits<double>::infinity() : diff[order[k-1]];
                const double hi=k==nf ? std::numeric_limits<double>::infinity() : diff[order[k]];
                if(candidate>=lo&&candidate<=hi){current=candidate;found=true;break;}
                if(k<nf){const int f=order[k];const double dr=dx/sp[f]-dx/sm[f];rr+=dr;rd+=dr*diff[f];}
            }
            if(!found) {
                // Rare rounded interval endpoints can leave no candidate when
                // the root coincides with a zero-field breakpoint. Invert the
                // same continuous monotone voltage function, without floors.
                const auto voltage_error=[&](double j){double v=j*rb-p[3];for(int f=0;f<nf;++f)v+=dx*(j-diff[f])/(j>=diff[f]?sp[f]:sm[f]);return v;};
                double low=std::min(0.0,*std::min_element(diff.begin(),diff.end()))-(std::abs(p[3])+1.0)/rb;
                double high=std::max(0.0,*std::max_element(diff.begin(),diff.end()))+(std::abs(p[3])+1.0)/rb;
                if(!(voltage_error(low)<=0.0&&voltage_error(high)>=0.0))throw std::runtime_error("Cannot bracket electrical current");
                double best=std::numeric_limits<double>::infinity();
                for(int it=0;it<120;++it){const double mid=.5*low+.5*high;const double residual=voltage_error(mid);
                    if(std::abs(residual)<best){current=mid;best=std::abs(residual);}
                    if(mid==low||mid==high||residual==0.0)break;
                    if(residual<0.0)low=mid;else high=mid;}
            }
            rr=rb;rd=0;for(int f=0;f<nf;++f){sf[f]=current>=diff[f]?sp[f]:sm[f];const double r=dx/sf[f];rr+=r;rd+=r*diff[f];}
            current=(p[3]+rd)/rr;
            for(int f=0;f<nf;++f)sf[f]=current>=diff[f]?sp[f]:sm[f];
            double volt=current*rb,maxout=0;
            std::fill(flux.begin(),flux.end(),0.0);
            for(int f=0;f<nf;++f) {
                const double E=(current-diff[f])/sf[f];volt+=E*dx;double actual=0;
                for(int s=0;s<5;++s){const double w=z[s]*FRT*D[s]*E;const int k=s*nf+f;
                    const double cl=c[s*N+f],cr=c[s*N+f+1];
                    if((cl==0.0&&recl[k]!=0.0)||(cr==0.0&&recr[k]!=0.0))throw std::runtime_error("Nonzero migration reconstruction from a zero donor");
                    const double fl=p[18]==1.0?1.0:(cl>0.0?recl[k]/cl:0.0);
                    const double fr=p[18]==1.0?1.0:(cr>0.0?recr[k]/cr:0.0);
                    left[k]=D[s]/dx+std::max(w,0.0)*fl;right[k]=D[s]/dx+std::max(-w,0.0)*fr;
                    const double J=left[k]*c[s*N+f]-right[k]*c[s*N+f+1];flux[s*(N+1)+f+1]=J;actual+=z[s]*J;}
                m[6]=std::max(m[6],std::abs(F*actual-current));
            }
            m[7]=std::max(m[7],std::abs(volt-p[3]));
            for(int s=0;s<5;++s)for(int i=0;i<N;++i){double out=0;if(i<nf)out+=left[s*nf+i];if(i>0)out+=right[s*nf+i-1];maxout=std::max(maxout,out);}
            const double dt=std::min(limit,p[9]*theta*dx/maxout);
            if(!(dt>0.0)||!std::isfinite(dt)||elapsed+dt==elapsed)throw std::runtime_error("Invalid or stalled CFL step");
            if(current>=0.0){flux[0]=current/F;flux[(N+1)+N]=-current/F;}
            else{flux[N+1]=-current/F;flux[N]=current/F;}
            for(int s=0;s<5;++s)for(int i=0;i<N;++i){const int k=s*N+i;next[k]=c[k]+dt/(theta*dx)*(flux[s*(N+1)+i]-flux[s*(N+1)+i+1]);
                if(next[k]<0.0||!std::isfinite(next[k]))throw std::runtime_error("Transport positivity failed");m[3]=std::min(m[3],next[k]);}
            return std::pair<double,double>(dt,current);
        };
        const auto react_state = [&](double* c,double* cs,double* cp) {
            for(int i=0;i<N;++i){const double total=c[2*N+i]+cs[i]+cp[i],charge=c[i]-c[N+i]+2.0*c[2*N+i];
                chem.solve(total,charge,c[i],c[i],c[N+i],c[2*N+i],cs[i],cp[i],m);}
        };
        const auto inspect_state = [&](const double* c,const double* cs,const double* cp) {
            double inv[3]={0,0,0};
            for(int i=0;i<N;++i) {
                double q=0,denom=0;for(int s=0;s<5;++s){const double value=c[s*N+i];if(value<0.0||!std::isfinite(value))throw std::runtime_error("Invalid reacted state");m[3]=std::min(m[3],value);q+=z[s]*value;denom+=std::abs(z[s]*value);}
                m[4]=std::max(m[4],std::abs(q));m[5]=std::max(m[5],std::abs(q)/denom);
                inv[0]+=c[2*N+i]+cs[i]+cp[i];inv[1]+=c[3*N+i];inv[2]+=c[4*N+i];
                const double h=c[i]/1000.0,oh=c[N+i]/1000.0,aq=c[2*N+i]/1000.0;
                const double base=chem.h50/h,ratio=p[7]==2.0?base*base:std::pow(base,p[7]);
                const double tot=c[2*N+i]+cs[i]+cp[i],ip=aq*oh*oh/p[5];
                m[11]=std::max(m[11],std::abs(h*oh/p[4]-1.0));
                m[12]=std::max(m[12],std::abs(cs[i]-ratio*c[2*N+i])/std::max(tot,1e-15));
                m[13]=std::max(m[13],std::max(ip-1.0,0.0));if(cp[i]>1e-12)m[13]=std::max(m[13],std::abs(ip-1.0));
            }
            for(int k=0;k<3;++k)m[8+k]=std::max(m[8+k],std::abs(inv[k]*scale/p[15+k]-1.0));
        };
        std::vector<double> stage1(5*N),sorb1(N),ppt1(N);
        while(elapsed<duration-1e-9) {
            double dt=0.0,faradaic_current=0.0,signed_current=0.0;
            if(p[20]==1.0) {
                const auto proposal=transport_proposal(c,std::min(p[8],duration-elapsed));
                dt=proposal.first;faradaic_current=std::abs(proposal.second);
                signed_current=proposal.second;
                m[22]=std::min(m[22],proposal.second);m[23]=std::max(m[23],proposal.second);
                std::copy(next.begin(),next.end(),c);
                react_state(c,cs,cp);
            } else {
                double limit=std::min(p[8],duration-elapsed);
                bool accepted=false;
                for(int attempt=0;attempt<=30;++attempt) {
                    const auto first=transport_proposal(c,limit);
                    dt=first.first;
                    std::copy(next.begin(),next.end(),stage1.begin());
                    std::copy(cs,cs+N,sorb1.begin());std::copy(cp,cp+N,ppt1.begin());
                    react_state(stage1.data(),sorb1.data(),ppt1.data());
                    inspect_state(stage1.data(),sorb1.data(),ppt1.data());
                    const auto second=transport_proposal(stage1.data(),dt);
                    if(second.first<dt*(1.0-1e-12)) {
                        // Restart from U0 using the stricter stage CFL. The
                        // trial states never overwrite the accepted solution.
                        limit=.99*second.first;m[21]+=1.0;continue;
                    }
                    for(int k=0;k<5*N;++k)c[k]=.5*(c[k]+next[k]);
                    for(int i=0;i<N;++i){cs[i]=.5*(cs[i]+sorb1[i]);cp[i]=.5*(cp[i]+ppt1[i]);}
                    react_state(c,cs,cp);
                    faradaic_current=.5*(std::abs(first.second)+std::abs(second.second));
                    signed_current=.5*(first.second+second.second);
                    m[22]=std::min({m[22],first.second,second.second});
                    m[23]=std::max({m[23],first.second,second.second});
                    accepted=true;break;
                }
                if(!accepted)throw std::runtime_error("RK2 stage CFL did not converge after 30 retries");
            }
            inspect_state(c,cs,cp);
            m[14]+=faradaic_current*dt/F;m[1]=std::min(m[1],dt);m[2]=std::max(m[2],dt);
            m[24]+=signed_current*dt;
            elapsed+=dt;m[0]+=1.0;m[15]=elapsed;
        }
        return 0;
    } catch(const std::exception& e) {
        if(err&&errlen>0)std::snprintf(err,errlen,"At chunk t=%.17g s: %s",elapsed,e.what());
        return 1;
    }
}

// Isolated equilibrium entry point permits direct accelerated/log validation.
extern "C" int equilibrate_native(int N,const double* total,const double* charge,
        const double* hint,const double* p,double* result,double* metrics,char* err,int errlen) {
    try {
        Chem chem(p);std::fill(metrics,metrics+21,0.0);
        for(int i=0;i<N;++i){
            if(!std::isfinite(total[i])||!std::isfinite(charge[i])||total[i]<0.0)
                throw std::runtime_error("Invalid equilibrium invariant");
            chem.solve(total[i],charge[i],hint[i],result[i],result[N+i],result[2*N+i],result[3*N+i],result[4*N+i],metrics);
        }
        return 0;
    } catch(const std::exception& e){if(err&&errlen>0)std::snprintf(err,errlen,"%s",e.what());return 1;}
}
