#!/bin/bash
for file in /usr/inputs/* ;do
    echo $file
    [ -f $file ] && oink --$1 "$file" -z $2 --no
done