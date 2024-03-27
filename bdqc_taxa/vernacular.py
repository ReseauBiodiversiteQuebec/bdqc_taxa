from . import gbif
from . import bryoquel
from . import cdpnq
from . import taxa_ref
from . import wikidata

# ACCEPTED_DATA_SOURCE = [
#     'Integrated Taxonomic Information System (ITIS)',
#     'Catalogue of Life Checklist']

ACCEPTED_LANGUAGE = ['fra', 'eng']

GBIF_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family',
                'genus', 'species', 'subspecies', 'variety', 'form']

def rank_order(rank):
    if not rank:
        return 9999
    try:
        return GBIF_RANKS.index(rank.lower())
    except ValueError:
        return 9999

def initcap_vernacular(name):
    capped_words = ['Amérique', 'America', 'Europe', 'Ungava', 'Alléghanys', 'Oregon', 'Virginie', 'Virginia', 'New York', 'Alaska', 'Pennsylvanie', 'Pennsylvania' , 'Canada', 'Inde', 'India', 'Islande', 'Égypte', 'Egypt', 'Pacifique', 'Pacific', 'Atlantique', 'Atlantic', 'Fraser', 'Est', 'Ouest', 'Nord', 'Alep', 'Anadyr', 'Eames', 'Allen', 'Anna', 'Uhler', 'Audubon']
    capped_words_lower = [word.lower() for word in capped_words]
    out = name[0].upper() + name[1:].lower()
    out_list = out.split(' ')
    for i, word in enumerate(out_list):
        if word.lower() in capped_words_lower:
            out_list[i] = out_list[i].capitalize()
    out = ' '.join(out_list)
    return out


class Vernacular:
    def __init__(self,
                 name: str = '',
                 source: str = '',
                 source_taxon_key: str = '',
                 language: str = '',
                 rank: str = None,
                 rank_order: int = 9999):
        self._name = name
        self.source = source
        self.source_taxon_key = source_taxon_key
        self.language = language.lower() if language else None
        self.rank = rank.lower() if rank else None
        self.rank_order = rank_order

    @property
    def name(self):
        return initcap_vernacular(self._name)

    @classmethod
    def from_gbif(cls, gbif_key: int, rank: str = None):
        out = []
        gbif_results = gbif.Species.get_vernacular_name(gbif_key)
        for result in gbif_results:
            if result['language'] not in ACCEPTED_LANGUAGE:
                continue
            vernacular = cls(
                    name = result['vernacularName'],
                    source = result['source'],
                    language = result['language'],
                    source_taxon_key = result['sourceTaxonKey'],
                    rank = rank,
                    rank_order = rank_order(rank)
                )
            out.append(vernacular)
        # dict comprehension trick to get only unique objects
        out = list({str(vars(o)): o for o in out}.values())
        return out
    
    @classmethod
    def from_gbif_match(cls, name: str = '', **match_kwargs):
        taxa = gbif.Species.match(name = name, **match_kwargs)

        try:
            return cls.from_gbif(taxa['usageKey'], rank=taxa['rank'])
        except KeyError:
            return []

    @classmethod
    def from_bryoquel_match(cls, name: str = ''):
        taxa = bryoquel.match_taxa(name)
        out = []
        if taxa and taxa['vernacular_name_fr']:
            out = [*out, cls(
                    name = taxa['vernacular_name_fr'],
                    source = 'Bryoquel',
                    language = 'fra',
                    source_taxon_key = taxa['id'],
                    rank = taxa['taxon_rank'],
                    rank_order = rank_order(taxa['taxon_rank'])
                )]
        if taxa and taxa['vernacular_name_en']:
            out = [*out, cls(
                    name = taxa['vernacular_name_en'],
                    source = 'Bryoquel',
                    language = 'eng',
                    source_taxon_key = taxa['id'],
                    rank = taxa['taxon_rank'],
                    rank_order = rank_order(taxa['taxon_rank'])
                )]
        return out

    @classmethod
    def from_cdpnq_match(cls, name: str = ''):
        taxas = cdpnq.match_taxa(name)
        if not taxas:
            return []
        out = []
        for taxa in taxas:
            out.append(cls(
                name = taxa['vernacular_fr'],
                source = 'CDPNQ',
                language = 'fra',
                source_taxon_key = taxa['name'],
                rank = taxa['rank'],
                rank_order = rank_order(taxa['rank'])
            ))
        
        return out

    @classmethod
    def from_wikidata_match(cls, name: str = '', rank: str = None):
        LANGUAGE_DICT = {
            'fr': 'fra',
            'en': 'eng'
        }
                
        out = []
        try:
            result = wikidata.search_entities(name)[0]
        except IndexError:
            return out
        
        entity = wikidata.get_entities(result['id'], languages=['fr', 'en'])

        # For each language, we only keep the first result that is not a scientific name
        for language in LANGUAGE_DICT.keys():
            try:
                name_dicts = [entity['labels'][language]]
            except KeyError:
                name_dicts = []

            try:
                name_dicts += entity['aliases'][language]
            except KeyError:
                pass

            if not name_dicts:
                continue

            for name_dict in name_dicts:
                vernacular = cls(
                    name = name_dict['value'],
                    source = 'Wikidata',
                    language = LANGUAGE_DICT[language],
                    source_taxon_key = result['id'],
                    rank = rank,
                    rank_order = rank_order(rank)
                )
                if vernacular.name.lower() != name.lower():
                    out.append(vernacular)
                    break

        return out

    @classmethod
    def from_match(cls, name: str, rank: str = None, **match_kwargs):
        out = cls.from_gbif_match(name, **match_kwargs)

        # Get the first result rank to use as a fallback
        if not rank and out:
            rank = GBIF_RANKS[out[0].rank_order]

        out = [*out, *cls.from_bryoquel_match(name)]
        out = [*out, *cls.from_cdpnq_match(name)]
        out = [*out, *cls.from_wikidata_match(name, rank = rank)]
        return out