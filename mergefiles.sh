FILES=`find -name *~$1`
for f in $FILES
do
    FILE=`echo $f | cut -d'~' -f 1`
    cat $FILE~$1 $FILE~$2 > $FILE
    git add $FILE
done
git merge --continue
