from . import gbif
from . import bryoquel
from inspect import signature

# ACCEPTED_DATA_SOURCE = [
#     'Integrated Taxonomic Information System (ITIS)',
#     'Catalogue of Life Checklist']

ACCEPTED_LANGUAGE = ['fra', 'eng']

def initcap(sentence):
    sentence = sentence.lower()
    return " ".join(
        [word[0].upper() + word[1:].lower() for word in sentence.split(" ")]
    )

class Vernacular:
    def __init__(self,
                 name: str = '',
                 source: str = '',
                 source_taxon_key: str = '',
                 language: str = ''):
        for param in signature(self.__init__).parameters:
            setattr(self, param, eval(param))
        self.language = language.lower()
        self.name = name.lower()
    
    @classmethod
    def from_gbif(cls, gbif_key: int):
        out = []
        gbif_results = gbif.Species.get_vernacular_name(gbif_key)
        for result in gbif_results:
            if result['language'] not in ACCEPTED_LANGUAGE:
                continue
            out.append(
                cls(
                    name = result['vernacularName'],
                    source = result['source'],
                    language = result['language'],
                    source_taxon_key = result['sourceTaxonKey']
                )
            )
        # dict comprehension trick to get only unique objects
        out = list({str(vars(o)): o for o in out}.values())
        return out
    
    @classmethod
    def from_gbif_match(cls, name: str = '', **match_kwargs):
        species = gbif.Species.match(name = name, **match_kwargs)
        try:
            return cls.from_gbif(species['usageKey'])
        except KeyError:
            return [None]

    @classmethod
    def from_bryoquel_match(cls, name: str = ''):
        species = bryoquel.match_species(name)
        try:
            return [cls(
                name = species['vernacular_name_fr'],
                source = 'Bryoquel',
                language = 'fra',
                source_taxon_key = species['id']
            ), cls(
                name = species['vernacular_name_en'],
                source = 'Bryoquel',
                language = 'eng',
                source_taxon_key = species['id']
            )]
        except KeyError:
            return [None]

    @classmethod
    def get(cls, name: str = None, gbif_key = None, **match_kwargs):
        if gbif_key:
            out = cls.from_gbif(gbif_key)
        elif name:
            out = cls.from_gbif_match(name = name, **match_kwargs)
        
        if name:
            out.extend(cls.from_bryoquel_match(name))
        
        return out