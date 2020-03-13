WHERE=$1

if [ -z "$1" ]
then
    export where=$PWD
else
    export where=$1
fi

rm $where/ran*
rm $where/*.log
rm $where/*.err
rm $where/*.out
rm $where/*_fort*


