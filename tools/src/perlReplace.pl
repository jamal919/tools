#!/usr/bin/perl

# usefull script to apply global multi-line changes

sub giveLine{
  ($package, $filename, $line)=caller;
  return $line;
  }

$configPath="../config.h.in";
$clip="";
$fileCount=0;

FILE: for $p (@ARGV){#for all files given as argument
  if($p eq "-h" or $p eq "--help"){
    print
      "USE\n".
      "  $0 file1 [file2 ...]\n".
      "\n".
      "OPTIONS\n".
      "  -h,--help  print this help and exit\n".
      "  -t  do not use cliboard for all following files (default)\n".
      "  -x  use X clipboard as custom expression to evaluate\n".
      "";
    exit 0;
    }
  
  if($p eq "-t"){
    $clip="";
    next;
    }
  
  if($p eq "-x"){
    $clip=`xclip -o`;
    print "Evaluating:\n$clip\n";
    next;
    }
  
  $fileCount++;
  
  if(not $p=~/\.(def|[ch](pp)?)$/ or $p=~/^(triangle.(h|cpp)|tricall.cpp|gshhs.h)$/){
    print "skipping $p\n";next
    }
  #read the file
  open($f,"<",$p);
  $_=join "",<$f>;
  $previous=$_;
  close($f);
  if(/\$Id: .*\$\s*?\n/){
    warn "skipping $p BECAUSE OF VERSION CONTROL SIGNATURE";next
    }
  
#------------------------------------------------------------------------------
  @tocheck=();
  @changed=();
  @maybe=();
  
#------------------------------------------------------------------------------
  # CHECKS
  if(0){
    #descriptions
    #file:///usr/share/doc/packages/doxygen/html/commands.html#cmdbrief says
    #A brief description ends when a blank line or another sectioning command is encountered.
    /print_help.*^\s*int main/sm or next;
    ($brief)=/[\\@]brief (.*?)\n\s*(?:\n|\\)/s;
    ($desc)=/".*DESCRIPTION.*".*\n.*"(?:\\n)? *(.*)"/;
    $brief or $desc or next;
    print "\n$p:\n$brief\n$desc\n";
    next;
    }
  
  # double include
  if(0){
    @_=sort /(#include [^ \n]+)/mg;
    for $i (1..$#_){
      if($_[$i-1] eq $_[$i]){
        warn "$p:$_[$i]";
        }
      }
    next;
    }
  
  # deletion and reinitialisation of pointers
  if(0){
    /delete.*\n?.*NULL/ or next;
    print "$p\n";
    next;
    }
  
#------------------------------------------------------------------------------
  # SPACES
  if(1){
    #delete redundant spaces between returned type and function name
    #can mess-up comments
    #s|(^ *(?:extern +)?[_a-z][_a-z0-9]+)(?<!do)(?<!else)  +(\*?[_a-z][_a-z0-9]* *\()|$1 $2|gmi and push @changed,giveLine;
    
    #ensure 2 empty lines between function bodies
    s|(\})(?: *\n){0,2}(/\*X+)|$1\n\n\n$2|gm and push @changed,giveLine;
    
    #delete redundant tabs
    $i=0;
    while(s|^((?:        )*) {0,7}\t|$1        |gm){$i++}
    $i and push @changed,giveLine."x$i";
    
    # delete redundant spaces
    # - between if() or for() and {}
    s/(\b(?:if|for)\b[^\n]+\))\s*( )\s+\{/$1$2\{/g and push @changed,giveLine;
    # - at the end of lines
    s|([^ \n\/])[ \t]+$|$1|mg and push @changed,giveLine;
    # - before template
    #s|^ +(template)|$1|gm and push @changed,giveLine;
    # - after '}'
    s|(\n( +)  }\n\2) +\n|$1\n|mg and push @changed,giveLine;
    
    #add spaces
    s|(XXX\*/\n) *([a-z])|$1\n$2|gi and push @changed,giveLine;
    s|\)(\n/\*XXX)|)\n$1|gi and push @maybe,giveLine;
    s|(\S)<<([^=\s])|$1 << $2|g and push @changed,giveLine;
    
    #remove spaces
    s|(#define[^\n]+\n)\n(/\*XXX)|$1$2|gi and push @maybe,giveLine;
    s| +(\n/\*X+\*/)|$1|gi and push @changed,giveLine;
    s|(/\*X+\*/ *\n) +\n|$1\n|gi and push @changed,giveLine;
    }
  
#------------------------------------------------------------------------------
  # COMMENTS
  $dash_x69="-"x69;
  $dash_x76="-"x76;
  $dash_x78="-"x78;
  $X_x76="X"x76;
  $star_x79="*"x79;
  if(1){
    #ensure template
    #/print_help.*^\s*int main/sm or next;s/\n\n(\s*?[\\@]brief)/\n\n<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->\n$1/ and push @changed,giveLine;
    s|^/\*-{1,75}\*/ *$|/*$dash_x76*/|mg and push @changed,giveLine;
    s|^/\*-{1,75} *$|/*$dash_x78|mg and push @changed,giveLine;
    s|^/\*X{1,75}\*/ *$|/*$X_x76*/|mgi and push @changed,giveLine;
    s|^/\*X{77,}\*/ *$|/*$X_x76*/|mgi and push @changed,giveLine;
    s|^/\*{10,78} *$|/$star_x79|mg and push @changed,giveLine;
    s|^\*{10,78}/ *$|$star_x79/|mg and push @changed,giveLine;
    
    s|\@\*/ *$|@ */|mg and push @changed,giveLine;
    
    #ensure proper titles
    if(1){
      s|^(?:/\*-+\*/\n)?( *)/// *<h1> *(?:it +)?([^\n]*?) *</h1> *\n|/*$dash_x69*//**<h1>\n$1$2 </h1>*/\n|gmi and push @changed,giveLine;
      }
    else{
      s|^( *)/// *<h1> *(?:it +)?([^\n]*?) *</h1> *\n|/*$dash_x76--\n$1$2 */\n|gmi and push @changed,giveLine;
      }
    }
  
#------------------------------------------------------------------------------
  if(1){
    #delete wrong include
    s|^ *#include <projects.h> *\n||mg and push @changed,giveLine;
    
    #change wrong include
    /#include <stdlib.h>/ and s|^ *#include <malloc.h> *\n||mg and push @changed,giveLine;
    s|#include <malloc.h>|#include <stdlib.h>|g and push @changed,giveLine;
    
    #delete redundant include
    /#include "spectrum.h"/ and s|^ *#include "tides.h" *\n||mg and push @changed,giveLine;
    /#include "poc-netcdf.hpp"/ and s|^ *#include "netcdf-proto.h" *\n||mg and push @changed,giveLine;
    
    #delete useless include
    #/(?:\bFILE|scanf|printf)\b/ or s|^ *#include <stdio.h> *\n||mg and push @changed,giveLine;
    /\b(?:malloc|free)\b/ or s|^ *#include <stdlib.h> *\n||mg and push @changed,giveLine;
    #/\bc(?:in|out|err|log)\b/ or s|^ *#include <iostream> *\n||mg and push @changed,giveLine;
    #/(?:filebuf|fstream)/ or s|^ *#include <fstream> *\n||mg and push @changed,giveLine;
    #/string(?:buf|stream)/ or s|^ *#include <sstream> *\n||mg and push @changed,giveLine;
    /\bvector\b *</ or s|^ *#include <vector> *\n||mg and push @changed,giveLine;
    /\blist\b *</ or s|^ *#include <list> *\n||mg and push @changed,giveLine;
    /[^\/]gsl_/ or s|^ *#include <gsl/gsl_.*\n||mg and push @changed,giveLine;
    /\b([fl]?stat|mkdir)\b *\(/ or s|^ *#include <sys/stat.h> *\n||mg and push @changed,giveLine;
    /\bu(int|short|long)\b */ or s|^ *#include <sys/types.h> *\n||mg and push @changed,giveLine;
    #kdevelop madness
    /\b(pj_|PJ\b)/ or s|^ *#include ["<](?:lib_proj4/)?lib_proj.h[">] *\n||mg and push @changed,giveLine;
    /\b(?:zapsc|zaparr)\b/ or s|^ *#include <zapper.h> *\n||mg and push @changed,giveLine;
    
    if($configDefsList eq "" and -r $configPath){
      open($f,"<",$configPath);
      $configDefsList=join "|",map {chomp;s/#undef // ? $_ : ()} <$f>;
      close($f);
      }
    /\b(?:$configDefsList)\b/ or s|^ *#include <config.h> *\n||mg and push @changed,giveLine;
    
    # delete redundant quotes in compilation warnings
    s|(#warning )"([^"]*)"|$1$2|g and push @changed,giveLine;
    
    # more lines on templates
    s|(^ *(?:template *<[^>]+>\s*)?([a-z][a-z_0-9<> *]*?)\s*\b[a-z][a-z_0-9]*\s*\([^)/;]+\)\s*/\*X+\*/\n\{)\s*return ([a-z][^;]+;)\s*\}|$1\n  $2 result;\n  \n  result=$3\n  \n  return result;\n}|msgi and push @changed,giveLine;
    }
  
  if(@names=(/(?:^|[^ 0-9A-Za-z_\.]|[^\/][\/\*])\s*\b(access|getopt|truncate|unlink|read|write|chdir)\b\s*\(/g)){#check include
    unless(/#include +<unistd.h>/){
      warn "$p MAY BADLY MISS #include <unistd.h> because of ".join(' ',@names);
      }
    }
  else{
    #s|^ *#include <unistd.h> *\n||mg and push @changed,giveLine;
    }
  
#------------------------------------------------------------------------------
  # UPDATES
  if(0){
    # put exitIfNull almost wherever necessary
    s/(\n *)if *\( *\( *([^=]+?) *= *(.*?) *\) *== *(?:NULL|0) *\) *\{[ \t]*\n *__ERR_BASE_LINE__\(""\);perror.*\n *w?exit.*\n *\}/$1exitIfNull($1  $2=$3$1  );/g and push @changed,giveLine;
    s/(\n *)(([^=]+?) *= *.*?) *; *\n* *if *\( *\3 *== *(?:NULL|0) *\) *\{[ \t]*\n *__ERR_BASE_LINE__\(""\);perror.*\n *w?exit.*\n *\}/$1exitIfNull($1  $2$1  );/g and push @changed,giveLine;
    }
  
  if($p =~ m|(.*/)functions.cpp$|){
    # replace exit by its wrapper
    s/\bexit\(/wexit(/g and push @changed,giveLine;
    }
  
  if(1){
    # replace geo_recale(,,180.) by its wrapper
    s/\bgeo_recale( *\(.+?), *180.0*f?\)/degree_recale$1)/g and push @changed,giveLine;
    
    # replace MIN and MAX
    s/\bMIN( *)\(/min$1(/g and push @tocheck,giveLine;
    s/\bMAX( *)\(/max$1(/g and push @tocheck,giveLine;
    s/( *)([^=]{6,}?) *= *(min|max) *\( *(?:\2 *,([^\)]+)|([^,]+), *\2)\) *;/$1update$3(&$2,$4$5);/g and push @tocheck,giveLine;
    
    # remove annoying and useless static_cast
    s/static_cast *< *double *> *\( *([+-]?[0-9\.]+(?:[eE][+-]?[0-9]+)?) *\)/$1/g and push @changed,giveLine;
    s/static_cast *< *float *> *\( *([+-]?[0-9\.]+(?:[eE][+-]?[0-9]+)?) *\)/$1f/g and push @changed,giveLine;
    s/static_cast *< *int *> *\( *([+-]?[0-9]+) *\)/$1/g and push @changed,giveLine;
    if(/static_cast/){
      warn "*** $p STILL HAS static_cast ***\n";
      }
    
    # remove annoying (<type>)
    s/\( *double *\) *([+-]?[0-9\.]+(?:[eE][+-]?[0-9]+)?)/$1/g and push @changed,giveLine;
    s/\( *float *\) *([+-]?[0-9\.]+(?:[eE][+-]?[0-9]+)?f?)/$1f/g and push @changed,giveLine;
    s/\( *const +char *\* *\) *"/"/g and push @changed,giveLine;
    s/^( *)\( *void *\) */$1/mg and push @tocheck,giveLine;
    
    if(0){
      # remove annoying e+00
      $r="(?:[^A-Za-z_][0-9]+(?:\.[0-9]*)?|[0-9]*\.[0-9]+)";
      s/($r[eE][+-]?)0+([0-9])/$1$2/g and push @changed,giveLine;
      s/($r[eE])\+/$1/g and push @changed,giveLine;
      s/($r)[eE][+-]?0/$1/g and push @changed,giveLine;
      }
    
    # and redundant non-standard d suffix
    # s/([+-]?[0-9]+(?:.[0-9]*)?(?:[eE][+-]?[0-9]+)?)d/$1/g and push @changed,giveLine; # BUGGY
    
    # update julian_day()
    s/julian_day\((.+?).month,\1.day,\1.year\)/julian_day($1)/mg and push @changed,giveLine;
    
    # days to second constant
    if(/#include +"(?:tides|spectrum).h"/){
      s|/3600(?:\.0*)?/24(?:\.0*)?|/d2s|g and push @changed,giveLine;
      s|\(?3600(?:\.0*)?\*24(?:\.0*)?\)?|d2s|g and push @changed,giveLine;
      }
    
    # replace wordy compare
    s/\.compare *\( *(".*?") *\) *== *0/==$1/g and push @changed,giveLine;
    
    # break wordy error traps
    s/^((?:\/\/)?[\t ]+)if\s*\(\s*\(\s*((\S+)\s*=\s*[_a-z][a-z_0-9]*\s*\((?:[^\)]*\(\))*[^\)]*\))\s*\)\s*(==|!=)\s*(?:NULL|0|NC_NOERR)\s*\)/$1$2;\n$1if($3$4\x30)/mgi and push @changed,giveLine;
    }
  
  if(1){
    # remove redundant {}
    s/\)\s*\{\s*((?:__)?(?:(?:STD)?(?:ERR|OUT)_BASE_LINE|TRAP_ERR_EXIT|NC_CHKERR_LINE_FILE)(?:__)?\([^;]+?\);)\s*\}/) $1/g and push @changed,giveLine;
    s/\{\s*((?:__)?(?:(?:STD)?(?:ERR|OUT)_BASE_LINE|TRAP_ERR_EXIT|NC_CHKERR_LINE_FILE)(?:__)?\([^;]+?\);)\s*\}/ $1/g and push @changed,giveLine;
    
    # use __func__, as per info:/gcc/Function%20Names
    s/__FUNCTION__/__func__/g and push @changed,giveLine;
    
    # use STDERR_BASE_LINE_FUNC
    s/__ERR_BASE_LINE__\("(?:%s:?)([^:].*?)",__func__/STDERR_BASE_LINE_FUNC("$1"/g and push @changed,giveLine;
    
    # err... no, STDERR_BASE_LINE_FUNC
    s/__(ERR|OUT)_BASE_LINE(_FUNC)?__/STD$1_BASE_LINE$2/g and push @changed,giveLine;
    s/__TRAP_ERR_(EXIT|RETURN)__/TRAP_ERR_$1/g and push @changed,giveLine;
    s/__(NC_[A-Z_]+|CHKERR_LINE_FILE)__/$1/g and push @changed,giveLine;
    
    # use TRAP_ERR_EXIT
    s/(^ *)status *= *([_a-z][a-z_0-9]*|(?:[+-] *)?[0-9]+?) *;\n\1(check_error\( *)status( *, *[^,]+, *__LINE__, *__FILE__, *1\);)/$1$3$2$4/gmi and push @changed,giveLine;
    s/(^ *)sprintf\(([_a-z][a-z_0-9]*),(.*?)\);\n\1check_error\((.+?), *\2, *__LINE__, *__FILE__, *1\);/$1TRAP_ERR_EXIT($4,$3);/gmi and push @changed,giveLine;
    s/(^ *)sprintf\(([_a-z][a-z_0-9]*),(.*?)\);\n\1check_error\((.+?), *\2, *__LINE__, *__FILE__, *0\);/$1STDERR_BASE_LINE($3);/gmi and push @changed,giveLine;
    s/\bcheck_error\( *([_a-zA-Z].+?),(.*?)(" *), *__LINE__, *__FILE__, *1\);/if($1!=0) TRAP_ERR_EXIT($1,$2\\n$3);/gm and push @changed,giveLine;
    s/\bcheck_error\((.*?)(" *), *__LINE__, *__FILE__, *1\);/TRAP_ERR_EXIT($1\\n$2);/g and push @changed,giveLine;
    # and NC_TRAP_ERROR
    s/\bnc_check_error\( *([_a-z][a-z_0-9]*), *__LINE__, *__FILE__, *1\);/NC_TRAP_ERROR(wexit,$1,1,"");/g and push @changed,giveLine;
    
    # clean wordy NC_TRAP_ERROR constructs
    s/(^ *)if\(status!=0\) *\{?\n *(NC_TRAP_ERROR\([^)]+?\);) *(\n +\})?/$1$2/gm and push @changed,giveLine;
    }
  
  s/\{\s*[_a-z][_a-z0-9]*=([_a-z][_a-z0-9]*);\s*if\(([_a-z][_a-z0-9]*)\)\{?\s*__NC_CHKERR_LINE_FILE__\(\1,(.+?)\);\s*\}?\s*(exit|return)\s*\(?\s*\1\s*\)?;\s*\}/__NC_TRAP_ERROR__($5,$1,$2,$3);/gi and push @tocheck,giveLine;
  s/\{\s*if\(([_a-z][_a-z0-9]*)\)\{?\s*__NC_CHKERR_LINE_FILE__\(([_a-z][_a-z0-9]*),(.+?)\);\s*\}?\s*(exit|return)\s*\(?\s*\2\s*\)?;\s*\}/__NC_TRAP_ERROR__($4,$2,$1,$3);/gi and push @tocheck,giveLine;
  s/\{\s*__NC_CHKERR_LINE_FILE__\(([_a-z][_a-z0-9]*),(.+?)\);\s*(exit|return)\s*\(?\s*\1\s*\)?;\s*\}/__NC_TRAP_ERROR__($3,$1,$2);/gi and push @changed,giveLine;
  
  #s|(int [^;]*\bstatus\b[^;=]*; *\n)(?: *\n)*( *return *\(?status\)?;)|$1  __CHKERR_LINE_FILE__(ENOEXEC,"not finished");exit(ENOEXEC);\n$2|g and push @changed,giveLine;#Disable some bad functions
  
  if(1){#give line and file of exit
    s|\bprintf(\s*\(.*?\)\s*;\s*w?exit\s*\(.*?\)\s*;)|STDOUT_BASE_LINE$1|g and push @changed,giveLine;
    s/\b(?:fprintf\s*\(\s*stderr\s*,|STDERR_BASE_LINE\()(.*?)\)\s*;\s*w?exit\s*\((.*?)\)\s*;/TRAP_ERR_EXIT($2,$1);/g and push @changed,giveLine;
    s|(\s+)(perror\s*\(.*?\)\s*;\s*w?exit\s*\(.*?\)\s*;)|$1STDERR_BASE_LINE("");$2|g and push @changed,giveLine;
    
    $changedbelow=0;
    sub replacement($$$$){my ($before,$spaces,$exitor,$exitval)=@_;
      $before=~/_LINE\b|_LINE_|print_help/ and
        return "$before$spaces$exitor";
      warn "************* CHANGING!!!";
      warn "<$before>";
      warn "<$spaces>";
      warn "<$exitor><$exitval>";
      $changedbelow=1;
      return $before.$spaces."TRAP_ERR_EXIT($exitval".',"exiting\n");';
      }
    s|([^\n]+)(\s*)(\bw?exit\s*\((.*?)\);)|replacement($1,$2,$3,$4)|ge;$changedbelow and push @changed,giveLine;
    
    s|\bSTDERR_BASE_LINE\((".*?")\);\s*w?exit\s*\((.*?)\)\s*;|TRAP_ERR_EXIT($2,$1);|g and push @changed,giveLine;
    }
  
#------------------------------------------------------------------------------
  # CUSTOM
  if($clip){
    $e=eval($clip); warn $@ if $@; $e and push @changed,giveLine;
    print "Clipboard expression returned $e\n";
    }
  
#------------------------------------------------------------------------------
  #s|See the main function for how this works <!-- A LINK TO main\(\) WILL NOT LINK TO THE MAIN OF THE RIGHT SOURCE ! -->\nand print_help\(\) for how to use this.|<!-- A LINK TO main\(\) or print_help\(\) WILL NOT LINK TO THE RIGHT SOURCE ! -->\nSee the main function for how this works\nand the print_help function for how to use this.| and push @changed,giveLine;#Doxygen madness
  #$deadDefine="DEFINE_global_astro_angles";s|^ *#if ($deadDefine).*?^\s*#endif //\1\n+||mgs and push @changed,giveLine;#delete dead code
  s|(^ *gettimeofday)\s*\(\s*([^\),]+)\s*,\s*NULL\s*\)|$1($2)|gs and push @changed,giveLine;#update code
  
#------------------------------------------------------------------------------
  if($#maybe>=0){
    if($previous eq $_){
      print "treated $p at l.".join(";",@maybe)." but nothing changed.\n";
      }
    else{
      push @changed,@maybe;
      }
    }
  $#changed==-1 and $#tocheck==-1 and next;#skip update if nothing to replace
  #backup
  rename $p,"$p~";
  #update the file
  open($f,">",$p);
  print {$f} $_;
  close($f);
  print "treated $p at l.".join(";",@changed).".";
  if($#tocheck>=0){
    print " MAKE SURE IT COMPILES BECAUSE OF CHANGES AT l.".join(",",@tocheck)." !!!";
    }
  print " Diff: $p".'{~,}'."\n";
  push(@treated,$p);
  }

$#treated>=0 or die "None of the $fileCount files treated.\n";

print join(" ",1+$#treated."/$fileCount treated :",@treated)."\n";
