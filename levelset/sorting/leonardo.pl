#!/usr/bin/perl
$max=32;$mm=$max-1;

open A,">leonardo.cc";

$b=1;$c=-1;
print A "const int smb[$max]={";
foreach (0..$mm) {
    print A "$b".($_==$mm?"};\n":",");
    $d=$b+$c+1;$c=$b;$b=$d;
}

$b=1;$c=-1;
print A "const int smc[$max]={";
foreach (0..$mm) {
    print A "$c".($_==$mm?"};\n":",");
    $d=$b+$c+1;$c=$b;$b=$d;
}

$b=1;$c=-1;
print A "const int smd[$max]={";
foreach (0..$mm) {
    $d=$b-$c;
    print A "$d".($_==$mm?"};\n":",");
    $d=$b+$c+1;$c=$b;$b=$d;
}
