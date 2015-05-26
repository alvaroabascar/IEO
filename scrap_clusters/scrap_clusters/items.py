# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# http://doc.scrapy.org/en/latest/topics/items.html

import scrapy


class GeneItem(scrapy.Item):
    # define the fields for your item here like:
    cluster = scrapy.Field()
    index = scrapy.Field()
    pos = scrapy.Field()
    geneId = scrapy.Field()
    geneDescription = scrapy.Field()
