    if(j<1) {
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=0;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            prtxy=(phi[ij+2*m+2]-phi[ij+2]-phi[ij+2*m]+phi[ij])*d2f;
            phix=xsp*(mf*mf*(-phi[ij]*2*me+phi[ij+1]*2*me+prbx*(2*e-1))
                +f*(2-f)*(-phi[ij+m]*2*me+phi[ij+m+1]*2*me+prtx*(2*e-1))
                -f*mf*(-plty*2*me+prty*2*me+prtxy*(2*e-1)));
            temp1=phi[ij]*me*me+phi[ij+1]*e*(2-e)-prbx*e*me;
            temp2=phi[ij+m]*me*me+phi[ij+m+1]*e*(2-e)-prtx*e*me;
            temp3=plty*me*me+prty*e*(2-e)-prtxy*e*me;
            phiy=ysp*(-2*mf*temp1+2*mf*temp2+(2*f-1)*temp3);
            return mf*mf*temp1+f*(2-f)*temp2-f*mf*temp3;
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=m-2;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            pltxy=(phi[ij+2*m+1]-phi[ij+1]-phi[ij+2*m-1]+phi[ij-1])*d2f;
            phix=xsp*(mf*mf*(plbx*(1-2*e)-phi[ij]*2*e+phi[ij+1]*2*e)
                +f*(2-f)*(pltx*(1-2*e)-phi[ij+m]*2*e+phi[ij+m+1]*2*e)
                -f*mf*(pltxy*(1-2*e)-plty*2*e+prty*2*e));
            temp1=plbx*e*me+phi[ij]*(1-e*e)+phi[ij+1]*e*e;
            temp2=pltx*e*me+phi[ij+m]*(1-e*e)+phi[ij+m+1]*e*e;
            temp3=pltxy*e*me+plty*(1-e*e)+prty*e*e;
            phiy=ysp*(-2*mf*temp1+2*mf*temp2+(2*f-1)*temp3);
            return mf*mf*temp1+f*(2-f)*temp2-f*mf*temp3;
        } else {
            e-=double(i);
            me=1-e;
            ij=i;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            pltxy=(phi[ij+2*m+1]-phi[ij+1]-phi[ij+2*m-1]+phi[ij-1])*d2f;
            prtxy=(phi[ij+2*m+2]-phi[ij+2]-phi[ij+2*m]+phi[ij])*d2f;
            phix=xsp*(mf*mf*(plbx*me*(1-3*e)-phi[ij]*6*e*me+phi[ij+1]*6*e*me+prbx*e*(3*e-2))
                +f*(2-f)*(pltx*me*(1-3*e)-phi[ij+m]*6*e*me+phi[ij+m+1]*6*e*me+prtx*e*(3*e-2))
                -f*mf*(pltxy*me*(1-3*e)-plty*6*e*me+prty*6*e*me+prtxy*e*(3*e-2)));
            temp1=plbx*e*me*me+phi[ij]*me*me*(1+2*e)+phi[ij+1]*e*e*(3-2*e)-prbx*e*e*me;
            temp2=pltx*e*me*me+phi[ij+m]*me*me*(1+2*e)+phi[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me;
            temp3=pltxy*e*me*me+plty*me*me*(1+2*e)+prty*e*e*(3-2*e)-prtxy*e*e*me;
            phiy=ysp*(-2*mf*temp1+2*mf*temp2+(2*f-1)*temp3);
            return mf*mf*temp1+f*(2-f)*temp2-f*mf*temp3;
        }
    } else if(j>=n-2) {
        f-=double(n-2);
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=mn-2*m;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            prbxy=(phi[ij+m+2]-phi[ij-m+2]-phi[ij+m]+phi[ij-m])*d2f;
            phix=xsp*(f*mf*(-plby*2*me+prby*2*me+prbxy*(2*e-1))
                +(1-f*f)*(-phi[ij]*2*me+phi[ij+1]*2*me+prbx*(2*e-1))
                +f*f*(-phi[ij+m]*2*me+phi[ij+m+1]*2*me+prtx*(2*e-1)));
            temp0=plby*me*me+prby*e*(2-e)-prbxy*e*me;
            temp1=phi[ij]*me*me+phi[ij+1]*e*(2-e)-prbx*e*me;
            temp2=phi[ij+m]*me*me+phi[ij+m+1]*e*(2-e)-prtx*e*me;
            phiy=ysp*((1-2*f)*temp0-2*f*temp1+2*f*temp2);
            return f*mf*temp0+(1-f*f)*temp1+f*f*temp2;
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=mn-m-2;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            plbxy=(phi[ij+m+1]-phi[ij-m+1]-phi[ij+m-1]+phi[ij-m-1])*d2f;
            phix=xsp*(f*mf*(plbxy*(1-2*e)-plby*2*e+prby*2*e)
                +(1-f*f)*(plbx*(1-2*e)-phi[ij]*2*e+phi[ij+1]*2*e)
                +f*f*(pltx*(1-2*e)-phi[ij+m]*2*e+phi[ij+m+1]*2*e));
            temp0=plbxy*e*me+plby*(1-e*e)+prby*e*e;
            temp1=plbx*e*me+phi[ij]*(1-e*e)+phi[ij+1]*e*e;
            temp2=pltx*e*me+phi[ij+m]*(1-e*e)+phi[ij+m+1]*e*e;
            phiy=ysp*((1-2*f)*temp0-2*f*temp1+2*f*temp2);
            return f*mf*temp0+(1-f*f)*temp1+f*f*temp2;
        } else {
            e-=double(i);
            me=1-e;
            ij=mn-2*m+i;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            plbxy=(phi[ij+m+1]-phi[ij-m+1]-phi[ij+m-1]+phi[ij-m-1])*d2f;
            prbxy=(phi[ij+m+2]-phi[ij-m+2]-phi[ij+m]+phi[ij-m])*d2f;
            phix=xsp*(f*mf*(plbxy*me*(1-3*e)-plby*6*e*me+prby*6*e*me+prbxy*e*(3*e-2))
                +(1-f*f)*(plbx*me*(1-3*e)-phi[ij]*6*e*me+phi[ij+1]*6*e*me+prbx*e*(3*e-2))
                +f*f*(pltx*me*(1-3*e)-phi[ij+m]*6*e*me+phi[ij+m+1]*6*e*me+prtx*e*(3*e-2)));
            temp0=plbxy*e*me*me+plby*me*me*(1+2*e)+prby*e*e*(3-2*e)-prbxy*e*e*me;
            temp1=plbx*e*me*me+phi[ij]*me*me*(1+2*e)+phi[ij+1]*e*e*(3-2*e)-prbx*e*e*me;
            temp2=pltx*e*me*me+phi[ij+m]*me*me*(1+2*e)+phi[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me;
            phiy=ysp*((1-2*f)*temp0-2*f*temp1+2*f*temp2);
            return f*mf*temp0+(1-f*f)*temp1+f*f*temp2;
        }
    } else {
        f-=double(j);
        mf=1-f;
        if(i<1) {
            me=1-e;
            ij=j*m;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            prbxy=(phi[ij+m+2]-phi[ij-m+2]-phi[ij+m]+phi[ij-m])*d2f;
            prtxy=(phi[ij+2*m+2]-phi[ij+2]-phi[ij+2*m]+phi[ij])*d2f;
            phix=xsp*(f*mf*mf*(-plby*2*me+prby*2*me+prbxy*(2*e-1))
                +mf*mf*(1+2*f)*(-phi[ij]*2*me+phi[ij+1]*2*me+prbx*(2*e-1))
                +f*f*(3-2*f)*(-phi[ij+m]*2*me+phi[ij+m+1]*2*me+prtx*(2*e-1))
                -f*f*mf*(-plty*2*me+prty*2*me+prtxy*(2*e-1)));
            temp0=plby*me*me+prby*e*(2-e)-prbxy*e*me;
            temp1=phi[ij]*me*me+phi[ij+1]*e*(2-e)-prbx*e*me;
            temp2=phi[ij+m]*me*me+phi[ij+m+1]*e*(2-e)-prtx*e*me;
            temp3=plty*me*me+prty*e*(2-e)-prtxy*e*me;
            phiy=ysp*(mf*(1-3*f)*temp0-6*f*mf*temp1+6*f*mf*temp2+f*(3*f-2)*temp3);
            return f*mf*mf*temp0+mf*mf*(1+2*f)*temp1+f*f*(3-2*f)*temp2-f*f*mf*temp3;
        } else if(i>=m-2) {
            e-=double(m-2);
            me=1-e;
            ij=(j+1)*m-2;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            plbxy=(phi[ij+m+1]-phi[ij-m+1]-phi[ij+m-1]+phi[ij-m-1])*d2f;
            pltxy=(phi[ij+2*m+1]-phi[ij+1]-phi[ij+2*m-1]+phi[ij-1])*d2f;
            phix=xsp*(f*mf*mf*(plbxy*(1-2*e)-plby*2*e+prby*2*e)
                +mf*mf*(1+2*f)*(plbx*(1-2*e)-phi[ij]*2*e+phi[ij+1]*2*e)
                +f*f*(3-2*f)*(pltx*(1-2*e)-phi[ij+m]*2*e+phi[ij+m+1]*2*e)
                -f*f*mf*(pltxy*(1-2*e)-plty*2*e+prty*2*e));
            temp0=plbxy*e*me+plby*(1-e*e)+prby*e*e;
            temp1=plbx*e*me+phi[ij]*(1-e*e)+phi[ij+1]*e*e;
            temp2=pltx*e*me+phi[ij+m]*(1-e*e)+phi[ij+m+1]*e*e;
            temp3=pltxy*e*me+plty*(1-e*e)+prty*e*e;
            phiy=ysp*(mf*(1-3*f)*temp0-6*f*mf*temp1+6*f*mf*temp2+f*(3*f-2)*temp3);
            return f*mf*mf*temp0+mf*mf*(1+2*f)*temp1+f*f*(3-2*f)*temp2-f*f*mf*temp3;
        } else {
            e-=double(i);
            me=1-e;
            ij=j*m+i;
            plbx=(phi[ij+1]-phi[ij-1])*dxf;
            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;
            prbx=(phi[ij+2]-phi[ij])*dxf;
            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;
            plby=(phi[ij+m]-phi[ij-m])*dyf;
            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;
            plty=(phi[ij+2*m]-phi[ij])*dyf;
            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;
            plbxy=(phi[ij+m+1]-phi[ij-m+1]-phi[ij+m-1]+phi[ij-m-1])*d2f;
            pltxy=(phi[ij+2*m+1]-phi[ij+1]-phi[ij+2*m-1]+phi[ij-1])*d2f;
            prbxy=(phi[ij+m+2]-phi[ij-m+2]-phi[ij+m]+phi[ij-m])*d2f;
            prtxy=(phi[ij+2*m+2]-phi[ij+2]-phi[ij+2*m]+phi[ij])*d2f;
            phix=xsp*(f*mf*mf*(plbxy*me*(1-3*e)-plby*6*e*me+prby*6*e*me+prbxy*e*(3*e-2))
                +mf*mf*(1+2*f)*(plbx*me*(1-3*e)-phi[ij]*6*e*me+phi[ij+1]*6*e*me+prbx*e*(3*e-2))
                +f*f*(3-2*f)*(pltx*me*(1-3*e)-phi[ij+m]*6*e*me+phi[ij+m+1]*6*e*me+prtx*e*(3*e-2))
                -f*f*mf*(pltxy*me*(1-3*e)-plty*6*e*me+prty*6*e*me+prtxy*e*(3*e-2)));
            temp0=plbxy*e*me*me+plby*me*me*(1+2*e)+prby*e*e*(3-2*e)-prbxy*e*e*me;
            temp1=plbx*e*me*me+phi[ij]*me*me*(1+2*e)+phi[ij+1]*e*e*(3-2*e)-prbx*e*e*me;
            temp2=pltx*e*me*me+phi[ij+m]*me*me*(1+2*e)+phi[ij+m+1]*e*e*(3-2*e)-prtx*e*e*me;
            temp3=pltxy*e*me*me+plty*me*me*(1+2*e)+prty*e*e*(3-2*e)-prtxy*e*e*me;
            phiy=ysp*(mf*(1-3*f)*temp0-6*f*mf*temp1+6*f*mf*temp2+f*(3*f-2)*temp3);
            return f*mf*mf*temp0+mf*mf*(1+2*f)*temp1+f*f*(3-2*f)*temp2-f*f*mf*temp3;
        }
    }
