# -*- coding: utf-8 -*-
import scrapy

from scrap_clusters.items import GeneItem

class ExampleSpider(scrapy.Spider):
    name = "clusters"
    allowed_domains = ["bcf.isb-sib.ch"]
    start_urls = [
        'http://bcf.isb-sib.ch/projects/cancer/glio/Results-56231395/rows/G%d_using.html' %i 
        for i in range(420)]

    def parse(self, response):
        clustername = response.url.split("/")[-1].replace("_using.html","")
        for sel in response.xpath('//table[2]/tr'):
                geneItem = GeneItem()
                geneItem['cluster'] = clustername
                geneItem['index'] = sel.xpath('td[1]/text()').extract()
                geneItem['pos'] = sel.xpath('td[2]/text()').extract()
                geneItem['geneId'] = sel.xpath('td[3]/text()').extract()
                geneItem['geneDescription'] = sel.xpath('td[5]/text()').extract()
                yield geneItem
