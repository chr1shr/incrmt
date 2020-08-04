#include "heap.hh"
#include "levelset.hh"

template<class d_option>
heap<d_option>::heap(levelset &ils,int &imem,int *&ihp) : ls(ils),
    m(ils.m),n(ils.n),mn(ils.mn),dx(ils.dx),dy(ils.dy),
    xsp(ils.xsp),ysp(ils.ysp),
    dx2(ils.dx2),dy2(ils.dy2),dxbydy(dx/dy),dx2bydy2(dxbydy*dxbydy),bpc(ils.bpc),
    mem(imem), hp(ihp), s(ils.s), bp(ils.bp), phi(ils.phi), w(1), z(ils,s) {}

template<class d_option>
void heap<d_option>::add_memory() {
    mem<<=1;printf("scale up to %d\n",mem);
    int *php(new int[mem]);
    for(int i=1;i<w;i++) php[i]=hp[i];
    delete [] hp;
    hp=php;
}

template<class d_option>
void heap<d_option>::stack_test(int i,int ij) {
    if(i>0) check_gridpoint(i-1,ij-1);
    if(i<m-1) check_gridpoint(i+1,ij+1);
    if(ij>=m) check_gridpoint(i,ij-m);
    if(ij<mn-m) check_gridpoint(i,ij+m);
}

template<class d_option>
void heap<d_option>::stack_test2(int i,int ij) {
    if(i>0) check_gridpoint2(i-1,ij-1);
    if(i<m-1) check_gridpoint2(i+1,ij+1);
    if(ij>=m) check_gridpoint2(i,ij-m);
    if(ij<mn-m) check_gridpoint2(i,ij+m);
}

template<class d_option>
void heap<d_option>::check_gridpoint(int i,int ij) {
    if(s[ij]==z.q) {
        s[ij]=z.p;
        add(ij,calc(i,ij));
    }
}

template<class d_option>
void heap<d_option>::check_gridpoint2(int i,int ij) {
    if(s[ij]==z.p&&bp[ij]!=bpc) {
        bp[ij]=bpc;
        if(w==mem) add_memory();
        hp[w++]=ij;
    }
}

template<class d_option>
void heap<d_option>::calc_tentative() {
    int c,i,ij;
    for(c=1;c<w;c++) {
        ij=hp[c];i=ij%m;
        phi[ij]=calc(i,ij);
    }
}

template<class d_option>
void heap<d_option>::add(int ij,double tphi) {
    if(w==mem) add_memory();
    phi[ij]=tphi;
    hp[w++]=ij;
}

template<class d_option>
void heap<d_option>::build_heap() {
    int ij,ij2,ij3;
    double tphi,tphi2;
    int tc,c,bc=(w-1)>>1;
    if((w&1)==0) {
        while(bc>=1) {
            tc=bc;
            c=tc<<1;
            while(c<w) {
                ij=hp[tc];tphi=phi[ij];
                ij2=hp[c];tphi2=phi[ij2];
                ij3=hp[c+1];
                if(z.cm(tphi,tphi2)) {
                    if (z.cm(tphi,phi[ij3])) break;else {c++;hp[tc]=ij3;}
                } else {
                    if (z.cm(phi[ij3],tphi2)) {c++;hp[tc]=ij3;}
                    else hp[tc]=ij2;
                }
                hp[c]=ij;tc=c;c=tc<<1;
            }
            bc--;
        }
    } else {
        int fw=w-1;
        while(bc>=1) {
            tc=bc;c=tc<<1;
            while(c<fw) {
                ij=hp[tc];tphi=phi[ij];
                ij2=hp[c];tphi2=phi[ij2];
                ij3=hp[c+1];
                if(z.cm(tphi,tphi2)) {
                    if (z.cm(tphi,phi[ij3])) break;else {c++;hp[tc]=ij3;}
                } else {
                    if (z.cm(phi[ij3],tphi2)) {c++;hp[tc]=ij3;}
                    else hp[tc]=ij2;
                }
                hp[c]=ij;tc=c;c=tc<<1;
            }
            if (c==fw) {
                ij=hp[tc];ij2=hp[c];
                if(z.cm(phi[ij2],phi[ij])) {hp[tc]=ij2;hp[c]=ij;}
            }
            bc--;
        }
    }
    for(c=1;c<w;c++) bp[hp[c]]=c;
}

template<class d_option>
void heap<d_option>::check_heap() {
    int bc,c=2;
    while(c<w) {
        bc=c>>1;
        if(z.cm(phi[hp[c]],phi[hp[bc]])) printf("Error from %d down to %d\n",c,bc);
        c++;
    }
}

template<class d_option>
void heap<d_option>::introduce(int ij,double tphi) {
    int c=w,bc=c>>1,ijt=hp[bc];
    while(bc>=1&&z.cm(tphi,phi[ijt])) {
        bp[ijt]=c;
        hp[c]=hp[bc];
        c=bc;bc=c>>1;
        ijt=hp[bc];
    }
    bp[ij]=c;
    phi[ij]=tphi;
    hp[c]=ij;
    w++;
}

template<class d_option>
void heap<d_option>::trickle(int c,int ij,double tphi) {
    int bc=c>>1,ijt=hp[bc];
    if(bc>=1&&z.cm(tphi,phi[ijt])) {
        do {
            bp[ijt]=c;
            hp[c]=hp[bc];
            c=bc;bc=c>>1;
            ijt=hp[bc];
        } while(bc>=1&&z.cm(tphi,phi[ijt]));
        bp[ij]=c;
        hp[c]=ij;
    }
    phi[ij]=tphi;
}

template<class d_option>
void heap<d_option>::update(int i,int ij) {
    if (s[ij]==z.q) {
        s[ij]=z.p;
        if(w==mem) add_memory();
        introduce(ij,calc(i,ij));
    } else if (s[ij]==z.p) {
        trickle(bp[ij],ij,calc(i,ij));
    }

}

template<class d_option>
void heap<d_option>::scan(int i,int ij) {
    if (s[ij]==z.p) trickle(bp[ij],ij,calc(i,ij));
}

template<class d_option>
double heap<d_option>::calc(int i,int ij) {
    double a,b,c,d;
    const double ti=2.0/3.0;
    if(ij>=m&&z.in(ij-m)) {
        if(ij<mn-m&&z.in(ij+m)) {
            if(i<m-1&&z.in(ij+1)) {
                if(i>0&&z.in(ij-1)) {
                    if(look_d(ij,a)) {
                        if(look_u(ij,c)) {
                            d=z.cm(a,c)?a:c;return threeway1(look_l(i,ij,a),look_r(i,ij,c),true,a,c,d);
                        } else {
                            if (look_l(i,ij,b)) {
                                if(look_r(i,ij,d)) {
                                    b=z.cm(b,d)?b:d;
                                    return threeway2(true,true,false,b,a,c);
                                } else return fourway(d,c,b,a);
                            } else {
                                if(look_r(i,ij,d)) return fourway(b,c,d,a);
                                else {
                                    b=z.cm(b,d)?b:d;
                                    return threeway2(false,true,false,b,a,c);
                                }
                            }
                        }
                    } else {
                        if(look_u(ij,c)) {
                            if (look_l(i,ij,b)) {
                                if(look_r(i,ij,d)) {
                                    b=z.cm(b,d)?b:d;
                                    return threeway2(true,true,false,b,c,a);
                                } else return fourway(d,a,b,c);
                            } else {
                                if(look_r(i,ij,d)) return fourway(b,a,d,c);
                                else {
                                    b=z.cm(b,d)?b:d;
                                    return threeway2(false,true,false,b,c,a);
                                }
                            }
                        } else {
                            d=z.cm(a,c)?a:c;return threeway1(look_l(i,ij,a),look_r(i,ij,c),false,a,c,d);
                        }
                    }
                } else {
                    return threeway2(look_r(i,ij,b),look_u(ij,a),look_d(ij,c),b,a,c);
                }
            } else {
                if(i>0&&z.in(ij-1)) return threeway2(look_l(i,ij,b),look_u(ij,a),look_d(ij,c),b,a,c);
                else {
                    c=look_d(ij,a)?z.ct(a,dy)*ti:z.ct(a,dy);
                    d=look_u(ij,b)?z.ct(b,dy)*ti:z.ct(b,dy);
                    return z.cm(c,d)?c:d;
                }
            }
        } else {
            if(i<m-1&&z.in(ij+1)) {
                if(i>0&&z.in(ij-1)) return threeway1(look_l(i,ij,a),look_r(i,ij,c),look_d(ij,b),a,c,b);
                else return corner(look_r(i,ij,a),look_d(ij,b),a,b);
            } else {
                if(i>0&&z.in(ij-1)) return corner(look_l(i,ij,a),look_d(ij,b),a,b);
                else return look_d(ij,a)?z.ct(a,dy)*ti:z.ct(a,dy);
            }
        }
    } else {
        if(ij<mn-m&&z.in(ij+m)) {
            if(i<m-1&&z.in(ij+1)) {
                if(i>0&&z.in(ij-1)) return threeway1(look_l(i,ij,a),look_r(i,ij,c),look_u(ij,b),a,c,b);
                else return corner(look_r(i,ij,a),look_u(ij,b),a,b);
            } else {
                if(i>0&&z.in(ij-1)) return corner(look_l(i,ij,a),look_u(ij,b),a,b);
                else return look_u(ij,a)?z.ct(a,dy)*ti:z.ct(a,dy);
            }
        } else {
            if(i<m-1&&z.in(ij+1)) {
                if(i>0&&z.in(ij-1)) {
                    c=look_l(i,ij,a)?z.ct(a,dx)*ti:z.ct(a,dx);
                    d=look_r(i,ij,b)?z.ct(b,dx)*ti:z.ct(b,dx);
                    return z.cm(c,d)?c:d;
                } else return look_r(i,ij,a)?z.ct(a,dx)*ti:z.ct(a,dx);
            } else {
                if(i>0&&z.in(ij-1)) return look_l(i,ij,a)?z.ct(a,dx)*ti:z.ct(a,dx);
                else return 0;
            }
        }
    }
}

template<class d_option>
double heap<d_option>::fourway(double &a,double &b,double &c,double &d) {
    double e=2*(c-a),f=2*(d-b),g,h;
    if (z.cm(e,f)) {
        g=e-a;h=e-b;
        if (z.cm(g,0)) g=0;
        if (z.cm(h,0)) h=0;
        if (g*g+h*h*dx2bydy2<dx2) {
            g=1.5*f-c;h=f-b;
            if (z.cm(g,0)) g=0;
            if (z.cm(h,0)) h=0;
            if (g*g+h*h*dx2bydy2<dx2) return corner(true,true,c,d);
            else return corner(true,false,c,b);
        } else return corner(false,false,a,b);
    } else {
        g=f-a;h=f-b;
        if (z.cm(g,0)) g=0;
        if (z.cm(h,0)) h=0;
        if (g*g+h*h*dx2bydy2<dx2) {
            g=e-a;h=1.5*e-d;
            if (z.cm(g,0)) g=0;
            if (z.cm(h,0)) h=0;
            if (g*g+h*h*dx2bydy2<dx2) return corner(true,true,c,d);
            else return corner(false,true,a,d);
        } else return corner(false,false,a,b);
    }
}

template<class d_option>
double heap<d_option>::threeway1(bool h1,bool h2,bool v,double &a,double &c,double &b) {
    double d;
    if (h1) {
        if (h2) {
            d=z.cm(a,c)?a:c;return corner(true,v,d,b);
        } else {
            d=2*(a-c);
            double e=d-c,f=(v?1.5*d:d)-b;
            if (z.cm(e,0)) e=0;
            if (z.cm(f,0)) f=0;
            return e*e+f*f*dx2bydy2<dx2?
                corner(true,v,a,b):
                corner(false,v,c,b);
        }
    } else {
        if (h2) {
            d=2*(c-a);
            double e=d-a,f=(v?1.5*d:d)-b;
            if (z.cm(e,0)) e=0;
            if (z.cm(f,0)) f=0;
            return e*e+f*f*dx2bydy2<dx2?
                corner(true,v,c,b):
                corner(false,v,a,b);
        } else {
            d=z.cm(a,c)?a:c;return corner(false,v,d,b);
        }
    }
}

template<class d_option>
double heap<d_option>::threeway2(bool h,bool v1,bool v2,double &b,double &a,double &c) {
    double d;
    if (v1) {
        if (v2) {
            d=z.cm(a,c)?a:c;return corner(h,true,b,d);
        } else {
            d=2*(a-c);
            double e=d-c,f=(h?1.5*d:d)-b;
            if (z.cm(e,0)) e=0;
            if (z.cm(f,0)) f=0;
            return f*f+e*e*dx2bydy2<dx2?
                corner(h,true,b,a):
                corner(h,false,b,c);
        }
    } else {
        if (v2) {
            d=2*(c-a);
            double e=d-a,f=(h?1.5*d:d)-b;
            if (z.cm(e,0)) e=0;
            if (z.cm(f,0)) f=0;
            return f*f+e*e*dx2bydy2<dx2?
                corner(h,true,b,c):
                corner(h,false,b,a);
        } else {
            d=z.cm(a,c)?a:c;return corner(h,false,b,d);
        }
    }
}

template<class d_option>
double heap<d_option>::corner(bool h,bool v,double &a,double &b) {
    const double f2a=2.25*(dy2+dx2),f2b=2.25*dy2+dx2,f2c=dy2+2.25*dx2,f2d=dy2+dx2;
    const double b1=2.25*dx2>dy2?1+2.25*dx2bydy2:1+dy2/(2.25*dx2);
    const double b2=dx2>dy2?1+dx2bydy2:1+dy2/dx2;
    const double b3=dx2>2.25*dy2?1+dx2bydy2/2.25:1+2.25*dy2/dx2;
    const double i2a=dy2/f2a,i2b=dy2/f2b,i2c=dy2/f2c,i2d=dy2/f2d;
    double f;
    if (h) {
        if (v) {
            f=1.5*(a-b);f*=f;
            return f*b2<f2a?z.ct(a*1.5+b*1.5*dx2bydy2,dxbydy*sqrt(f2a-f))*i2a:corner_bail(true,true,a,b);
        } else {
            f=a-1.5*b;f*=f;
            return f*b3<f2b?z.ct(a*1.5+b*dx2bydy2,dxbydy*sqrt(f2b-f))*i2b:corner_bail(true,false,a,b);
        }
    } else {
        if (v) {
            f=1.5*a-b;f*=f;
            return f*b1<f2c?z.ct(a+b*1.5*dx2bydy2,dxbydy*sqrt(f2c-f))*i2c:corner_bail(false,true,a,b);
        } else {
            f=a-b;f*=f;
            return f*b2<f2d?z.ct(a+b*dx2bydy2,dxbydy*sqrt(f2d-f))*i2d:corner_bail(false,false,a,b);
        }
    }
}

template<class d_option>
double heap<d_option>::corner_bail(bool h,bool v,double &a,double &b) {
#ifdef CORNER_BAIL_WARNING
    printf("corner bail with %s %s %g %g\n",h?"2nd":"1st",v?"2nd":"1st",a,b);
#endif
    const double ti=2.0/3.0;
    double ht=h?z.ct(a,dx)*ti:z.ct(a,dx);
    double vt=v?z.ct(b,dy)*ti:z.ct(b,dy);
    return z.cm(ht,vt)?ht:vt;
}

template<class d_option>
bool heap<d_option>::look_l(int i,int ij,double &e) {
#ifndef FORCE_FIRST_ORDER
    if (i>1&&z.in(ij-2)) {
        if(z.cm(phi[ij-2],phi[ij-1])) {
            e=2*phi[ij-1]-0.5*phi[ij-2];return true;
        }
    }
#endif
    e=phi[ij-1];return false;
}

template<class d_option>
bool heap<d_option>::look_r(int i,int ij,double &e) {
#ifndef FORCE_FIRST_ORDER
    if (i<m-2&&z.in(ij+2)) {
        if(z.cm(phi[ij+2],phi[ij+1])) {
            e=2*phi[ij+1]-0.5*phi[ij+2];return true;
        }
    }
#endif
    e=phi[ij+1];return false;
}

template<class d_option>
bool heap<d_option>::look_d(int ij,double &e) {
#ifndef FORCE_FIRST_ORDER
    if (ij>=2*m&&z.in(ij-2*m)) {
        if(z.cm(phi[ij-2*m],phi[ij-m])) {
            e=2*phi[ij-m]-0.5*phi[ij-2*m];return true;
        }
    }
#endif
    e=phi[ij-m];return false;
}

template<class d_option>
bool heap<d_option>::look_u(int ij,double &e) {
#ifndef FORCE_FIRST_ORDER
    if (ij<mn-2*m&&z.in(ij+2*m)) {
        if(z.cm(phi[ij+2*m],phi[ij+m])) {
            e=2*phi[ij+m]-0.5*phi[ij+2*m];return true;
        }
    }
#endif
    e=phi[ij+m];return false;
}

template<class d_option>
void heap<d_option>::march(double limit) {
    int c,bc;
    int i,ij,ijt,ij2,ij3;
    double tphi,tphi2;
    while(w>1&&z.cm(phi[hp[1]],limit)) {
        ij=hp[1];
        z.add(ij);
        i=ij%m;;s[ij]=9;
        s[ij]=z.o;w--;
        ijt=hp[w];
        tphi=phi[ijt];
        bc=1;c=bc<<1;
        while(c+1<w) {
            ij2=hp[c];tphi2=phi[ij2];ij3=hp[c+1];
            if(z.cm(tphi2,tphi)) {
                if(z.cm(phi[ij3],tphi2)) {hp[bc]=ij3;bp[ij3]=bc;bc=c+1;}
                else {hp[bc]=ij2;bp[ij2]=bc;bc=c;}
            } else {
                if(z.cm(phi[ij3],tphi)) {hp[bc]=ij3;bp[ij3]=bc;bc=c+1;}
                else break;
            }
            c=bc<<1;
        }
        if (c+1==w) {
            ij2=hp[c];
            if(z.cm(phi[ij2],tphi)) {hp[bc]=ij2;bp[ij2]=bc;bc=c;}
        }
        hp[bc]=ijt;
        bp[ijt]=bc;
        if (ij>=2*m) {scan(i,ij-2*m);update(i,ij-m);} else if (ij>=m) update(i,ij-m);
        if (i>1) {scan(i-2,ij-2);update(i-1,ij-1);} else if (i>0) update(i-1,ij-1);
        if (i<m-2) {update(i+1,ij+1);scan(i+2,ij+2);} else if (i<m-1) update(i+1,ij+1);
        if (ij<mn-2*m) {update(i,ij+m);scan(i,ij+2*m);} else if (ij<mn-m) update(i,ij+m);
    }
}

template<class d_option>
void heap<d_option>::full_march() {
    int c,bc;
    int i,ij,ijt,ij2,ij3;
    double tphi,tphi2;
    while(w>1) {
        ij=hp[1];
        z.add(ij);
        i=ij%m;s[ij]=z.o;w--;
        ijt=hp[w];
        tphi=phi[ijt];
        bc=1;c=bc<<1;
        while(c+1<w) {
            ij2=hp[c];tphi2=phi[ij2];ij3=hp[c+1];
            if(z.cm(tphi2,tphi)) {
                if(z.cm(phi[ij3],tphi2)) {hp[bc]=ij3;bp[ij3]=bc;bc=c+1;}
                else {hp[bc]=ij2;bp[ij2]=bc;bc=c;}
            } else {
                if(z.cm(phi[ij3],tphi)) {hp[bc]=ij3;bp[ij3]=bc;bc=c+1;}
                else break;
            }
            c=bc<<1;
        }
        if (c+1==w) {
            ij2=hp[c];
            if(z.cm(phi[ij2],tphi)) {hp[bc]=ij2;bp[ij2]=bc;bc=c;}
        }
        hp[bc]=ijt;
        bp[ijt]=bc;
        if (ij>=2*m) {scan(i,ij-2*m);update(i,ij-m);} else if (ij>=m) update(i,ij-m);
        if (i>1) {scan(i-2,ij-2);update(i-1,ij-1);} else if (i>0) update(i-1,ij-1);
        if (i<m-2) {update(i+1,ij+1);scan(i+2,ij+2);} else if (i<m-1) update(i+1,ij+1);
        if (ij<mn-2*m) {update(i,ij+m);scan(i,ij+2*m);} else if (ij<mn-m) update(i,ij+m);
    }
}

template<class d_option>
double heap<d_option>::compute_eikonal(int i,int ij,int &stat) {
    double phid=0,phil=0,phir=0,phiu=0;
    if(ij>=m&&z.in(ij-m)) phid=(look_d(ij,phid)?1.5*phi[ij]-phid:phi[ij]-phid)*ysp;
    if(i>0&&z.in(ij-1)) phil=(look_l(i,ij,phil)?1.5*phi[ij]-phil:phi[ij]-phil)*xsp;
    if(i<m-1&&z.in(ij+1)) phir=(look_r(i,ij,phir)?1.5*phi[ij]-phir:phi[ij]-phir)*xsp;
    if(ij<mn-m&&z.in(ij+m)) phiu=(look_u(ij,phiu)?1.5*phi[ij]-phiu:phi[ij]-phiu)*ysp;
    if(z.cm(phiu,phid)) {
        if(z.cm(phid,0)) {stat=0;phid=0;}
        else stat=8;
    } else {
        if(z.cm(phiu,0)) {stat=0;phid=0;}
        else {stat=4;phid=phiu;}
    }
    if(z.cm(phir,phil)) {
        if(z.cm(phil,0)) phil=0;
        else stat|=2;
    } else {
        if(z.cm(phir,0)) phil=0;
        else {stat|=1;phil=phir;}
    }
    return sqrt(phil*phil+phid*phid);
}

dir_positive::dir_positive(levelset &ls_,int *&s) : ls(ls_),o(5),p(6),q(7),ss(s) {}
dir_negative::dir_negative(levelset &ls_,int *&s) : ls(ls_),o(2),p(1),q(0),ss(s) {}

void dir_positive::add(int ij) {ls.add_positive(ij);}
void dir_negative::add(int ij) {ls.add_negative(ij);}

// Explicit instantiation
template class heap<dir_positive>;
template class heap<dir_negative>;
