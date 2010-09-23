#!/usr/bin/perl

$tmp = (<STDIN>);
$tmp = " " . $tmp;
($trash, $num) = split(/\s+/, $tmp);
print $num; print "  ";
print "1"; print "\n";

while (<STDIN>) {
  chomp;
  $tmp = " " . $_;
  ($trash, $num, $v1, $v2, $v3) = split(/\s+/, $tmp);
  if ($v1 ne "Generated") {
    print $v1; print "  ";
    print $v2; print "  ";
    print $v3; print "\n";
  }
}
