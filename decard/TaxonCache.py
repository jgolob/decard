#!/usr/bin/env python
from Bio import Entrez
import json
import urllib
import re

class TaxonCache():

    def __init__(self, source='NCBI', url ="", email=""):
        self.__source__ = source
        self.__url__ = url
        self.__taxon_cache__ = {}
        
        self.rankOrders = {
            '0':  'kingdom',
            '1':  'phylum',
            '2':  'class',
            '3':  'order',
            '4':  'family',
            '5':  'genus',
            '6':  'species',
            'k':  'kingdom',
            'p':  'phylum',
            'c':  'class',
            'o':  'order',
            'f':  'family',
            'g':  'genus',
            's':  'species',
        }
        
        self.re_taxonomy_choices = [re.compile('D_(?P<rankID>\d+)__(?P<Name>.+)'),
                       re.compile('(?P<rankID>\w)__(?P<Name>.+)')
                       ]
    
        self.re_badNames = [
            re.compile('Incertae Sedis', re.IGNORECASE),
            re.compile('uncultured', re.IGNORECASE),
            re.compile('family', re.IGNORECASE),
            re.compile('unclassified', re.IGNORECASE),
            re.compile('unidentified', re.IGNORECASE),
            re.compile('human gut', re.IGNORECASE),
            re.compile('metagenome', re.IGNORECASE),
            re.compile('symbiont', re.IGNORECASE),
            re.compile('delta proteobacterium', re.IGNORECASE),
            re.compile('bacterium NLAE', re.IGNORECASE),
            re.compile('bacterium enrichment', re.IGNORECASE),
            re.compile('anaerobic bacterium', re.IGNORECASE),
            re.compile('\s+sp$', re.IGNORECASE),
            re.compile('\s+sp\.$', re.IGNORECASE),
            re.compile('\s+spp$', re.IGNORECASE),
            re.compile('\s+spp\.$', re.IGNORECASE),
            re.compile('^bacterium$', re.IGNORECASE),
            re.compile('\s+bacterium$', re.IGNORECASE),
            re.compile('^bacterium\s+', re.IGNORECASE),
            re.compile('drinking water', re.IGNORECASE),
            re.compile('phage', re.IGNORECASE),
            re.compile('denitrifying', re.IGNORECASE),
            re.compile('^iron$', re.IGNORECASE),
            re.compile('^bromate$', re.IGNORECASE),
            
        ]
        
        if email:
            Entrez.email = email
        
        if source=='custom':
            self.__lookup_id_method__ = self.__lookup_tax_id_custom__
            self.__lookup_name_method__ = self.__lookup_tax_name_custom__
        else:
            self.__lookup_id_method__ = self.__lookup_tax_id_NCBI__
            self.__lookup_name_method__ = self.__lookup_tax_name_NCBI__
        

        
    def __lookup_tax_id_custom__(self,tax_id):
        tax_id = int(tax_id)
        tax_h = urllib.urlopen(self.__url__+"id/"+str(tax_id)+"/")
        try:
            self.__taxon_cache__[tax_id]=json.load(tax_h)
            return self.__taxon_cache__[tax_id]
        except:
            return self.__lookup_tax_id_NCBI__(tax_id)
    
    def __lookup_tax_id_NCBI__(self,tax_id):
        tax_id = int(tax_id)
        entrez_h = Entrez.efetch(db='taxonomy', id = str(tax_id))
        
        try:
            self.__taxon_cache__[tax_id]=Entrez.read(entrez_h)[0]
            taxon = self.__taxon_cache__[tax_id]
            return taxon
        except:
            return None
    
    def lookup_tax_id(self, tax_id):
        # Convert to an int first
        tax_id = int(tax_id)
        # Try to look it up in our cache
        try:
            return self.__taxon_cache__[tax_id]
        except KeyError: # Not in our cache
            # look it up.
            return self.__lookup_id_method__(tax_id)
    
    # ------        
    def __lookup_tax_name_custom__(self,name):
        tax_h = urllib.urlopen(self.__url__+"name/"+name+"/")
        try:
            tax = json.load(tax_h)
            self.__taxon_cache__[tax['TaxId']]=tax
        except:
            return self.__lookup_tax_name_NCBI__(name)
        return tax
    
    def __lookup_tax_name_NCBI__(self,name):
        search_h = Entrez.esearch(db='taxonomy', term=name)
        search_result = Entrez.read(search_h)
        search_h.close()
    
        if int(search_result['Count']) > 0:
            return self.lookup_tax_id(int(search_result['IdList'][0]))
        else:
            return None
            
    def lookup_taxonomy(self, taxonomyList):
        """
            Input is a taxonomy in LIST format
        """
        tax = None
        if isinstance(taxonomyList,list):
            taxonomy_list = list(taxonomyList)
            
            while (tax == None and taxonomy_list):
                candidate_name = taxonomy_list.pop()
                if candidate_name:
                    tax = self.lookup_tax_name(candidate_name)
                
        return tax
    
    def lookup_tax_name(self, name):
        for tax_id in self.__taxon_cache__:
            if name.lower() == self.__taxon_cache__[tax_id]['ScientificName'].lower():
                return self.__taxon_cache__[tax_id]
        # Didn't find this name in our cache, next use our repo
        return self.__lookup_name_method__(name)
    # ----
    
    def is_bad_name(self, candidate_name):
        return len([r for r in self.re_badNames if r.search(candidate_name)]) > 0
        
    def get_best_name(self, taxonomyList):
        """
            This function is given a list of taxonomy [deepest,shallowest] and attempts to get the best sensible name, by stripping out things like 'uncultured bacterium'
            as found in the self.re_badNames list. Returns the first distinct name possible.
            Input: A taxononomy, with each rank broken out as a string in a list. Eg.
            ['Bacteria',
             'Firmicutes',
             'Bacilli',
             'Bacillales',
             'Family XII',
             'Exiguobacterium',
             'Exiguobacterium sp']
             
             Output: A string, with the best possible name. In this case 'Exiguobacterium sp' (or Exiguobacterium)
        """
        # Copy the input list, as we'll be popping from it.
        taxonomy_list = list(taxonomyList)
        best_name = ''
                
        while taxonomy_list and (best_name =='' or self.is_bad_name(best_name)):
            best_name = taxonomy_list.pop()
            best_name = best_name.replace("_"," ")
    
        return best_name

