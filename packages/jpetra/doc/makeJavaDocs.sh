#!/bin/sh

cd ../src
javadoc -d ../doc/html -linksource -author -private Jpetra Jpetra.MatrixMarketIO Jpetra.CcjSupport
