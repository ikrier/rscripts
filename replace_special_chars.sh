#!/bin/bash
#To replace weird hex ascii that annovar introduces in files :
sed -ri 's/\\x3b/\;/g;s/\\x3e/=/g' $1
