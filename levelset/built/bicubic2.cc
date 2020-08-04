    if(j<1) {
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=0;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            prtxy=(fe[ij+2*m+2]-fe[ij+2]-fe[ij+2*m]+fe[ij])*d2f;
            return mf*mf*(fe[ij]*me*me+fe[ij+1]*e*(2-e)-prbx*e*me)
                +f*(2-f)*(fe[ij+m]*me*me+fe[ij+m+1]*e*(2-e)-prtx*e*me)
                -f*mf*(plty*me*me+prty*e*(2-e)-prtxy*e*me);
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=m-2;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            pltxy=(fe[ij+2*m+1]-fe[ij+1]-fe[ij+2*m-1]+fe[ij-1])*d2f;
            return mf*mf*(plbx*e*me+fe[ij]*(1-e*e)+fe[ij+1]*e*e)
                +f*(2-f)*(pltx*e*me+fe[ij+m]*(1-e*e)+fe[ij+m+1]*e*e)
                -f*mf*(pltxy*e*me+plty*(1-e*e)+prty*e*e);
        } else {
            e-=double(i);
            me=1-e;
            ij=i;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            pltxy=(fe[ij+2*m+1]-fe[ij+1]-fe[ij+2*m-1]+fe[ij-1])*d2f;
            prtxy=(fe[ij+2*m+2]-fe[ij+2]-fe[ij+2*m]+fe[ij])*d2f;
            return mf*mf*(plbx*e*me*me+fe[ij]*me*me*(1+2*e)+fe[ij+1]*e*e*(3-2*e)-prbx*e*e*me)
                +f*(2-f)*(pltx*e*me*me+fe[ij+m]*me*me*(1+2*e)+fe[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me)
                -f*mf*(pltxy*e*me*me+plty*me*me*(1+2*e)+prty*e*e*(3-2*e)-prtxy*e*e*me);
        }
    } else if(j>=n-2) {
        f-=double(n-2);
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=mn-2*m;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            prbxy=(fe[ij+m+2]-fe[ij-m+2]-fe[ij+m]+fe[ij-m])*d2f;
            return f*mf*(plby*me*me+prby*e*(2-e)-prbxy*e*me)
                +(1-f*f)*(fe[ij]*me*me+fe[ij+1]*e*(2-e)-prbx*e*me)
                +f*f*(fe[ij+m]*me*me+fe[ij+m+1]*e*(2-e)-prtx*e*me);
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=mn-m-2;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            plbxy=(fe[ij+m+1]-fe[ij-m+1]-fe[ij+m-1]+fe[ij-m-1])*d2f;
            return f*mf*(plbxy*e*me+plby*(1-e*e)+prby*e*e)
                +(1-f*f)*(plbx*e*me+fe[ij]*(1-e*e)+fe[ij+1]*e*e)
                +f*f*(pltx*e*me+fe[ij+m]*(1-e*e)+fe[ij+m+1]*e*e);
        } else {
            e-=double(i);
            me=1-e;
            ij=mn-2*m+i;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            plbxy=(fe[ij+m+1]-fe[ij-m+1]-fe[ij+m-1]+fe[ij-m-1])*d2f;
            prbxy=(fe[ij+m+2]-fe[ij-m+2]-fe[ij+m]+fe[ij-m])*d2f;
            return f*mf*(plbxy*e*me*me+plby*me*me*(1+2*e)+prby*e*e*(3-2*e)-prbxy*e*e*me)
                +(1-f*f)*(plbx*e*me*me+fe[ij]*me*me*(1+2*e)+fe[ij+1]*e*e*(3-2*e)-prbx*e*e*me)
                +f*f*(pltx*e*me*me+fe[ij+m]*me*me*(1+2*e)+fe[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me);
        }
    } else {
        f-=double(j);
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=j*m;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            prbxy=(fe[ij+m+2]-fe[ij-m+2]-fe[ij+m]+fe[ij-m])*d2f;
            prtxy=(fe[ij+2*m+2]-fe[ij+2]-fe[ij+2*m]+fe[ij])*d2f;
            return f*mf*mf*(plby*me*me+prby*e*(2-e)-prbxy*e*me)
                +mf*mf*(1+2*f)*(fe[ij]*me*me+fe[ij+1]*e*(2-e)-prbx*e*me)
                +f*f*(3-2*f)*(fe[ij+m]*me*me+fe[ij+m+1]*e*(2-e)-prtx*e*me)
                -f*f*mf*(plty*me*me+prty*e*(2-e)-prtxy*e*me);
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=(j+1)*m-2;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            plbxy=(fe[ij+m+1]-fe[ij-m+1]-fe[ij+m-1]+fe[ij-m-1])*d2f;
            pltxy=(fe[ij+2*m+1]-fe[ij+1]-fe[ij+2*m-1]+fe[ij-1])*d2f;
            return f*mf*mf*(plbxy*e*me+plby*(1-e*e)+prby*e*e)
                +mf*mf*(1+2*f)*(plbx*e*me+fe[ij]*(1-e*e)+fe[ij+1]*e*e)
                +f*f*(3-2*f)*(pltx*e*me+fe[ij+m]*(1-e*e)+fe[ij+m+1]*e*e)
                -f*f*mf*(pltxy*e*me+plty*(1-e*e)+prty*e*e);
        } else {
            e-=double(i);
            me=1-e;
            ij=j*m+i;
            plbx=(fe[ij+1]-fe[ij-1])*dxf;
            pltx=(fe[ij+m+1]-fe[ij+m-1])*dxf;
            prbx=(fe[ij+2]-fe[ij])*dxf;
            prtx=(fe[ij+m+2]-fe[ij+m])*dxf;
            plby=(fe[ij+m]-fe[ij-m])*dyf;
            prby=(fe[ij+m+1]-fe[ij-m+1])*dyf;
            plty=(fe[ij+2*m]-fe[ij])*dyf;
            prty=(fe[ij+2*m+1]-fe[ij+1])*dyf;
            plbxy=(fe[ij+m+1]-fe[ij-m+1]-fe[ij+m-1]+fe[ij-m-1])*d2f;
            pltxy=(fe[ij+2*m+1]-fe[ij+1]-fe[ij+2*m-1]+fe[ij-1])*d2f;
            prbxy=(fe[ij+m+2]-fe[ij-m+2]-fe[ij+m]+fe[ij-m])*d2f;
            prtxy=(fe[ij+2*m+2]-fe[ij+2]-fe[ij+2*m]+fe[ij])*d2f;
            return f*mf*mf*(plbxy*e*me*me+plby*me*me*(1+2*e)+prby*e*e*(3-2*e)-prbxy*e*e*me)
                +mf*mf*(1+2*f)*(plbx*e*me*me+fe[ij]*me*me*(1+2*e)+fe[ij+1]*e*e*(3-2*e)-prbx*e*e*me)
                +f*f*(3-2*f)*(pltx*e*me*me+fe[ij+m]*me*me*(1+2*e)+fe[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me)
                -f*f*mf*(pltxy*e*me*me+plty*me*me*(1+2*e)+prty*e*e*(3-2*e)-prtxy*e*e*me);
        }
    }
