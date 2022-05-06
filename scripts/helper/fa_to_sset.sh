#BAZE=$(basename $1 .fa)
cat $1 | sed '/^>/d'  
