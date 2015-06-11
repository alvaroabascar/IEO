#!/bin/bash

rm clusters.csv
scrapy crawl clusters -o clusters.csv
