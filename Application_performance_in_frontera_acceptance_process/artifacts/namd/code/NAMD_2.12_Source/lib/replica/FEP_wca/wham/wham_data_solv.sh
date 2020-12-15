# To_change: The number of lines of date you will use
nLineUse=200000
echo "" > all-dat.tmp
tail -n $nLineUse ../output_solv/*/*.fepout >> all-dat.tmp

sed 's/ BREAK /\n/g' all-dat.tmp > org-dat.tmp
mv org-dat.tmp all-dat.tmp

head -n 1 ../output_solv/*/*.fepout | grep "^FEP_WCA_REP" | awk '{print $3}' | sort | uniq > all-cut2.list

grep "^FEP_ELEC" all-dat.tmp | awk '{print $2 "  " $3}' > chg.wham
grep "^FEP_WCA_DISP" all-dat.tmp | awk '{print $2 "  " $3}' > disp.wham


if [ ! -e "./extract_repul" ]; then
  g++ -O2 -o extract_repul extract-repul.cpp
fi

./extract_repul all-cut2.list all-dat.tmp 

rm all-dat.tmp


