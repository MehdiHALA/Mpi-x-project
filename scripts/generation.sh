for var in 10 25 50 75 100 200 300 500 750 1000
do
    ../build/stencil-basic $var $var | grep -oP 'GFLOPs.{0,13}' | cut -d= -f 2
done

