#!/bin/sh

cd ../src
javadoc -d ../doc/html-javadoc -linksource -author -private Jpetra Jpetra.MatrixMarketIO Jpetra.CcjSupport
