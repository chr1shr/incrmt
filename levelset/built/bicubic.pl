open B,">bicubic.cc";

@func=(["","mx*mx","x*(2-x)","-x*mx"],
    ["x*mx","(1-x*x)","x*x",""],
    ["x*mx*mx","mx*mx*(1+2*x)","x*x*(3-2*x)","-x*x*mx"]);

@dfunc=(["","-2*mx","2*mx","(2*x-1)"],
    ["(1-2*x)","-2*x","2*x",""],
    ["mx*(1-3*x)","-6*x*mx","6*x*mx","x*(3*x-2)"]);

@expr=(["plbxy","plby","prby","prbxy"],["plbx","phi[ij]","phi[ij+1]","prbx"],
    ["pltx","phi[ij+m]","phi[ij+m+1]","prtx"],["pltxy","plty","prty","prtxy"]);

foreach $k (0..2) {
    foreach $l (0..3) {
        $_=$func[$k][$l];
        $sgn[$k][$l]=m/^-/?"-":"+";
        $nsgn[$k][$l]=m/^-/?"-":"";
        s/^-//g;
        s/x/e/g;$e[$k][$l]=$_;
        s/e/f/g;$f[$k][$l]=$_;

        $_=$dfunc[$k][$l];
        $dsgn[$k][$l]=m/^-/?"-":"+";
        $ndsgn[$k][$l]=m/^-/?"-":"";
        s/^-//g;
        s/x/e/g;$de[$k][$l]=$_;
        s/e/f/g;$df[$k][$l]=$_;
    }
}

foreach $b (0,1,2) {

    if ($b==0) {
        print B "    if(j<1) {\n";
    } elsif ($b==1) {
        print B "    } else if(j>=n-2) {\n";
        print B "        f-=double(n-2);\n";
    } else {
        print B "    } else {\n";
        print B "        f-=double(j);\n";
    }
    print B "        mf=1-f;\n";

    foreach $a (0,1,2) {

        if ($a==0) {
            print B "        if(i<1) {\n";
        } elsif($a==1) {
            print B "        } else if(i>=m-2) {\n";
            print B "            e-=double(m-2);\n";
        } else {
            print B "        } else {\n";
            print B "            e-=double(i);\n";
        }
        print B "            me=1-e;\n";

        $c=3*$b+$a;
        if ($c==0) {$ij="0";}
        elsif ($c==1) {$ij="m-2";}
        elsif ($c==2) {$ij="i";}
        elsif ($c==3) {$ij="mn-2*m";}
        elsif ($c==4) {$ij="mn-m-2";}
        elsif ($c==5) {$ij="mn-2*m+i";}
        elsif ($c==6) {$ij="j*m";}
        elsif ($c==7) {$ij="(j+1)*m-2";}
        elsif ($c==8) {$ij="j*m+i";}

        print B "            ij=$ij;\n";

        if($a!=0) {
            print B "            plbx=(phi[ij+1]-phi[ij-1])*dxf;\n";
            print B "            pltx=(phi[ij+m+1]-phi[ij+m-1])*dxf;\n";
        }
        if($a!=1) {
            print B "            prbx=(phi[ij+2]-phi[ij])*dxf;\n";
            print B "            prtx=(phi[ij+m+2]-phi[ij+m])*dxf;\n";
        }
        if($b!=0) {
            print B "            plby=(phi[ij+m]-phi[ij-m])*dyf;\n";
            print B "            prby=(phi[ij+m+1]-phi[ij-m+1])*dyf;\n";
        }
        if($b!=1) {
            print B "            plty=(phi[ij+2*m]-phi[ij])*dyf;\n";
            print B "            prty=(phi[ij+2*m+1]-phi[ij+1])*dyf;\n";
        }
        if($b!=0&&$a!=0) {
            print B "            plbxy=(phi[ij+m+1]-phi[ij-m+1]-phi[ij+m-1]+phi[ij-m-1])*d2f;\n";
        }
        if($b!=1&&$a!=0) {
            print B "            pltxy=(phi[ij+2*m+1]-phi[ij+1]-phi[ij+2*m-1]+phi[ij-1])*d2f;\n";
        }
        if($b!=0&&$a!=1) {
            print B "            prbxy=(phi[ij+m+2]-phi[ij-m+2]-phi[ij+m]+phi[ij-m])*d2f;\n";
        }
        if($b!=1&&$a!=1) {
            print B "            prtxy=(phi[ij+2*m+2]-phi[ij+2]-phi[ij+2*m]+phi[ij])*d2f;\n";
        }
        $st=0;
        foreach $l (0..3) {
            next if $f[$b][$l] eq "";

            if ($st==0) {
                print B "            phix=xsp*($nsgn[$b][$l]$f[$b][$l]*(";$st=1;
            } else {
                print B ")\n                $sgn[$b][$l]$f[$b][$l]*(";
            }

            $st2=0;
            foreach $k (0..3) {
                next if $e[$a][$k] eq "";

                if($st2==0) {
                    print B "$ndsgn[$a][$k]";$st2=1;
                } else {
                    print B "$dsgn[$a][$k]";
                }

                print B "$expr[$l][$k]*$de[$a][$k]";
            }
        }
        print B "));\n";

        foreach $l (0..3) {
            next if $f[$b][$l] eq "";

            print B "            temp$l=";

            $st2=0;
            foreach $k (0..3) {
                next if $e[$a][$k] eq "";

                if($st2==0) {
                    print B "$nsgn[$a][$k]";$st2=1;
                } else {
                    print B "$sgn[$a][$k]";
                }

                print B "$expr[$l][$k]*$e[$a][$k]";
            }
            print B ";\n";
        }

        $st=0;
        foreach $l (0..3) {
            next if $f[$b][$l] eq "";

            if ($st==0) {
                print B "            phiy=ysp*($ndsgn[$b][$l]";$st=1;
            } else {
                print B "$dsgn[$b][$l]";
            }
            print B "$df[$b][$l]*temp$l";
        }
        print B ");\n";

        $st=0;
        foreach $l (0..3) {
            next if $f[$b][$l] eq "";

            if ($st==0) {
                print B "            return $nsgn[$b][$l]";$st=1;
            } else {
                print B "$sgn[$b][$l]";
            }
            print B "$f[$b][$l]*temp$l";
        }
        print B ";\n";
    }
    print B "        }\n";
}
print B "    }\n";
