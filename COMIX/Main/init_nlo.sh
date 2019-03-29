#!/bin/bash

dbtodir(){
  mkdir $1; fl=$(sqlite3 $1.db "select file from path");
  printf "Extracting $1.db ("$(echo $fl | wc -w)") ";
  for i in $fl; do printf "."; sqlite3 $1.db \
    "select content from path where file='$i'" > $1/$i;
  done; rm $1.db; echo " done"
}
dirtodb(){
  pf=$(find $1 -type f);
  printf "Compressing $1 ("$(echo $pf | wc -w)") ";
  ( echo -e "begin transaction;\ncreate table path(file,content);"
    for i in $pf; do fn=$(echo $i | sed 's|'$1'||g;s|^/||g');
    sed -e "s|'|''|g" -e "$ s|$|');|1" \
        -e "1 s|^|insert into path values('$fn','|1" $i;
    printf "." 1>&2; done; echo "commit;" ) | sqlite3 $1.db;
  rm -r $1; echo " done"
}

if test $# -lt 2; then 
  echo "usage: $0 <sherpa exe> <input file> [<process path>]";
  exit 1;
fi;
pd=$PWD; test -z "$3" || pd=$PWD/$3; echo $0": proc dir is '"$pd"'";
nt=$(mktemp -d -p$PWD); echo $0": temp dir is '"$nt"'";
tp=$(grep NLO_QCD $2 | sed -e's/.*NLO_QCD_Part[ \t]*\(\w*\).*/\1/g');
if test $tp = RS; then
  if ! grep -q "MEH_RSADD[ \t=]*0" $2; then
    echo "Input file must contain 'MEH_RSADD 0;'";
    rm -rf $nt; exit 1;
  fi;
fi;
sed -e's/}(run)/  INIT_ONLY 1;\n}(run)/g' < $2 > $2.$tp;
sed -e'/NLO_QCD/ d' < $2.$tp > $2.B;
test -z "$4" && cp -r $pd/Process/ $nt/ 2>&1;
$1 -f$2.B SHERPA_CPP_PATH=$nt;
dbtodir $nt/Process/Comix;
for i in $nt/Process/Comix/*.map; do
  if echo $i | grep -q 'QCD('$tp')'; then continue; fi
  if grep -q x $i; then
    sed -e's/ /__QCD('$tp') /g' $i > $i.tmp;
    mv $i.tmp $(echo $i | sed -e's/\(__NQ_.*\)[.]map/.map\1/g' \
      -e's/.map/__QCD('$tp').map/g' -e's/.map\(.*\)/\1.map/g');
  else
    test -f $(echo $i | sed -e's/.map/__QCD('$tp').map/g') && continue;
    sed -e'1 s/ /__QCD('$tp') /g' -e'1 s/$/__QCD('$tp')/g' $i > $i.tmp;
    if awk '{ if ($1!=$2) exit 1; exit 0; }' < $i; then rm $i.tmp;
    else mv $i.tmp $(echo $i | sed -e's/.map/__QCD('$tp').map/g'); fi;
  fi
done;
dirtodb $nt/Process/Comix;
$1 -f$2.$tp SHERPA_CPP_PATH=$nt;
echo -n $0": copying files ...";
cp -ur $nt/Process $pd;
echo " done";
rm -rf $nt;
