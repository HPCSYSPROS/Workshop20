num_repu=`ls -l repu_*.wham | wc -l`

if [ ! -e "./sort-data-wham" ]; then
  g++ -O2 -o sort-data-wham sort-data-fe.cpp
fi

if [ ! -e "./wham_fep" ]; then
  g++ -O2 -o wham_fep wham-fep.cpp
fi

for (( i=0; i<$num_repu; i++ ))
do
   nline=`wc -l repu_$i.wham | awk '{print $1}'`
   ./sort-data-wham repu_$i.wham repu_$i.wham $nline
   ./wham_fep repu_$i.wham repu_${i}_wham_fe > log-repu-$i.txt
done

nline=`wc -l chg.wham | awk '{print $1}'`
./sort-data-wham chg.wham chg.wham $nline
./wham_fep chg.wham chg_wham_fe > log-chg.txt

nline=`wc -l disp.wham | awk '{print $1}'`
./sort-data-wham disp.wham disp.wham $nline
./wham_fep disp.wham disp_wham_fe > log-disp.txt

