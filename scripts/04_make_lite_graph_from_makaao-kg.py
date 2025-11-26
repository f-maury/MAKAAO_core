# KG simplifier — notebook cell
from rdflib import Graph, Namespace, URIRef, RDF, RDFS, OWL
from rdflib.namespace import SKOS
from rdflib.term import Literal

# ---- Paths ----
IN_PATH  = "../kg/makg-core_v1.rdf"
OUT_PATH = "../kg/makg-core_v1_simplified.rdf"

# ---- Namespaces ----
MAK   = Namespace("http://makaao.inria.fr/kg/")
OBO   = Namespace("http://purl.obolibrary.org/obo/")
UNI   = Namespace("http://purl.uniprot.org/core/")
SIO   = Namespace("http://semanticscience.org/resource/")
BIOL  = Namespace("https://w3id.org/biolink/vocab/")
BAO   = Namespace("http://www.bioassayontology.org/bao#")
PROV  = Namespace("http://www.w3.org/ns/prov#")
ORDO  = Namespace("http://www.orpha.net/ORDO/")

# ---- Important URIs ----
CHEBI_23367   = URIRef(OBO["CHEBI_23367"]) # chebi : molecular entity
PROTEIN       = URIRef(UNI["Protein"])
TARGET        = URIRef(MAK["Target"])
AUTOIMMUNE_DZ = URIRef(MAK["AutoimmuneDisease"])
ORPHANET_C001 = URIRef(ORDO["Orphanet_C001"]) # orphanet : disease
EXCLUDE_CLASSES = {URIRef(MAK["Document"]), URIRef(MAK["Relation"])}

# ---- Whitelists ----
KEEP_PRED = {RDFS.label, SKOS.prefLabel, RDFS.subClassOf}
REL_PRED  = {
    SIO["SIO_001403"], # is_associated_with
    BIOL["biomarker_for"],
    BIOL["has_biomarker"],
    BAO["BAO_0000598"], # is_target_for
    BAO["BAO_0000211"], # has_target
    SIO["SIO_001279"], # has_phenotype
    SIO["SIO_001280"] # is_phenotype_of
}

# ---- Helpers ----
def tail(u: URIRef):
    return str(u).rsplit("/", 1)[-1] # keep onlu slug of URI

def is_uri(x):
    return isinstance(x, URIRef) # check if something is an URI

def is_prov(x):
    return is_uri(x) and str(x).startswith(str(PROV)) # check if something has a PROV: URI

def is_instance_uri(u: URIRef): # check if URI is an instance (Document_/Relation_ or ends with _instance) (check this works in every case ??)
    s = str(u)
    if s.endswith("_instance"):
        return True
    if s.startswith(str(MAK)):
        t = tail(u)
        if t.startswith("document_"):
            return True
        if t.startswith("r") and t[1:].isdigit():
            return True
    return False

def strip_instance(u: URIRef): # remove _instance at end of URI
    s = str(u)
    return URIRef(s[:-9]) if s.endswith("_instance") else None

def copy_labels_from_instance(src_inst: URIRef, dst_class: URIRef, outg: Graph, gsrc: Graph):
    for o in gsrc.objects(src_inst, RDFS.label):
        outg.add((dst_class, RDFS.label, o)) # copy labels of an instance of full KG, to a class in simplified KG
    for o in gsrc.objects(src_inst, SKOS.prefLabel): # copy labels of an instance of full KG, to a class in simplified KG
        outg.add((dst_class, SKOS.prefLabel, o))

def has_label(gx: Graph, node: URIRef): # check if an item has at least 1 label rdfs:label or skos:prefLabel
    for _ in gx.objects(node, RDFS.label):
        return True
    for _ in gx.objects(node, SKOS.prefLabel):
        return True
    return False

def cui_class_of(u: URIRef): # get the CUI URI, from a CUI instance in full KG
    s = str(u)
    if "/concept/C" in s:
        return URIRef(s[:-9]) if s.endswith("_instance") else URIRef(s)
    return None

# ---- Load source graph ----
g = Graph()
g.parse(IN_PATH) # load input KG

# Precompute types
types_by_node = {}
for s, o in g.subject_objects(RDF.type):
    if is_prov(s) or is_prov(o):
        continue
    types_by_node.setdefault(s, set()).add(o) # filter out PROV type items

# Promotions: local UniProt instances UP_*_instance -> UP_* (reuse labels)
promotion_UP_local = {
    n: strip_instance(n) for n in types_by_node # associate Uniprot instances to their classes URI
    if is_uri(n) and str(n).startswith(str(MAK)) and tail(n).startswith("UP_") and strip_instance(n) is not None # check if we can strip _instance
}
promotion_map = dict(promotion_UP_local)


def to_class(node: URIRef):
    """Resolve instance→class with explicit rules."""
    if not is_uri(node):
        return None
    if not is_instance_uri(node):
        return node
    t = tail(node)
    s = str(node)

    if node in promotion_map: #uniprot nistance to class dict
        return promotion_map[node]
    if t.startswith("CHEBI_"):
        chebi_id = t[:-9] if t.endswith("_instance") else t # remove _instance at end of chebi id
        return URIRef(str(OBO) + chebi_id) # reconstruct class chebi URI
    if "/concept/C" in s and s.endswith("_instance"): # if we have a CUI instance, we remove instance: we obtain a UMLS URI
        base = strip_instance(node)
        if base is not None:
            return base
    if t.startswith("orpha_") and t.endswith("_instance"): # if we detect an orphanet instance
        num = t[len("orpha_"):-len("_instance")] # we extract the orphanet number
        if num.isdigit():
            return URIRef(str(ORDO) + f"Orphanet_{num}") # we reconstruct the orphanet class URI with the number
    if t.startswith("positivity_") and t.endswith("_instance"): # if we detect a positivity instance
        ts = types_by_node.get(node, set())
        for u in ts: # for each type of that positivity instance
            if str(u).startswith(str(OBO)) and tail(u).startswith("HP_"): # we check if the type is an HP class, if yes we get that URI
                return u                  # e.g., HP_0034117
        return strip_instance(node) # else we just strip _instance
    if s.startswith(str(MAK)) and t.endswith("_instance"): # if we have a local MAK (like an AAb class) instance
        return strip_instance(node)

    # 7) No rule matched
    return None

def all_nodes(gx: Graph):
    nodes = set()
    for s, p, o in gx:
        nodes.add(s); nodes.add(o) #  get all subject-object pairs in a graph
    return nodes

# ---- Build simplified graph ----
out = Graph() # start new empty graph

def involves_prov(s,p,o): # check if any of the subject, predicate or object is a PROV URI
    return is_prov(s) or is_prov(p) or is_prov(o)

# 1) Copy whitelisted class-level triples, excluding Document/Relation
for s, p, o in g:
    if involves_prov(s, p, o): # we don't add to the new graph triples with PROV:
        continue
    if p not in KEEP_PRED: # we don't add in the new graph triples with non whitelisted predicates
        continue
    if isinstance(s, URIRef) and not is_instance_uri(s) and s not in EXCLUDE_CLASSES: # if a subject URI is not an instance, and not in excluded class: we add the triple
        out.add((s, p, o))

# Keep explicit owl:Class declarations, excluding banned classes
for s, o in g.subject_objects(RDF.type): # we check triple with rdf:type predicate only
    if is_prov(s) or is_prov(o): # skip it  if subject or object is PROV
        continue
    if o == OWL.Class and isinstance(s, URIRef) and not is_instance_uri(s) and s not in EXCLUDE_CLASSES: #  if object is Class, subject is not instance, not excluded class
        out.add((s, RDF.type, OWL.Class)) # we add it

# 2) Rewire selected predicates from instance→class
for s, p, o in g: # for each triple in source graph
    if p not in REL_PRED or involves_prov(s, p, o): # if triple not using an allowed predicate, or involves PROV, skip
        continue
    s_cls = to_class(s) # get class URI of subject
    o_cls = to_class(o) # get class URI of object
    if s_cls is None or o_cls is None: # if we can't get class URI, we skip
        continue
    if is_instance_uri(s_cls) or is_instance_uri(o_cls): # if we still have an instnce URI, we skip
        continue
    if s_cls in EXCLUDE_CLASSES or o_cls in EXCLUDE_CLASSES: # if either subject or object class is in excluded classes, we skip
        continue
    out.add((s_cls, p, o_cls)) #else, we add the rewired triple to the new graph

# 2b) Reuse labels for UniProt promotions on their classes
for inst, cls in promotion_map.items():
    copy_labels_from_instance(inst, cls, out, g) # copy labels from instance to Uniprot class in new graph

# ---- Build CUI classification sets from source graph (direct + reified) ----
def build_cui_sets(gsrc: Graph):
    cui_target, cui_disease, cui_all = set(), set(), set()

    def note(node):
        c = cui_class_of(node) # get CUI class from instance, add it to set of all CUI classes
        if c is not None:
            cui_all.add(c)
        return c

    # direct triples,we have to check using predicate, othewise we would not know if that CUI is a Target or a Disease
    for s, p, o in gsrc: # for all triples in orginal graph
        if p in (BAO["BAO_0000211"], BAO["BAO_0000598"]): # if predicate is "is_target_of" or "has_target"
            for n in (s, o):
                c = note(n) # we try to get CUI class from subject or object
                if c is not None:
                    cui_target.add(c)
        elif p == SIO["SIO_001403"]: # if predicate is "is_related_to"
            for n in (s, o):
                c = note(n)
                if c is not None: # we try to get CUI class from subject or object
                    cui_disease.add(c) # if we got a CUI class, it is necessarily form at disease (the other end of that relation would ba an aab)


    return cui_all, cui_target, cui_disease

cui_all, cui_target, cui_disease = build_cui_sets(g) # call function on original graphn to get CUI classes that are DIeases, and those that are Targets

# 3) Enforce constraints
def force_constraints(outg: Graph):
    nodes = all_nodes(outg) # get all subject-object pairs in new graph

        # CHEBI normalization
    chebis = {n for n in nodes if isinstance(n, URIRef) and tail(n).startswith("CHEBI_")} # we get all triples involving chebi classes

    # Non-root CHEBI_* ⊑ CHEBI_23367 only
    for c in chebis - {CHEBI_23367}: # for all these chebi class, except root chebi
        for _, _, o in list(outg.triples((c, RDFS.subClassOf, None))): # we check their superclass
            if o != CHEBI_23367:
                outg.remove((c, RDFS.subClassOf, o)) # and we rmeove them, if not root chebi
        outg.add((c, RDFS.subClassOf, CHEBI_23367)) # we make them subclasses of root chebi only
        outg.add((c, RDF.type, OWL.Class))

    # Root: CHEBI_23367 ⊑ Target only
    for _, _, o in list(outg.triples((CHEBI_23367, RDFS.subClassOf, None))): # we check superclasses of root chebi
        if o != TARGET:
            outg.remove((CHEBI_23367, RDFS.subClassOf, o)) # we remove them if not target
    outg.add((CHEBI_23367, RDFS.subClassOf, TARGET)) # we add Target as only superclass of root chebi
    outg.add((CHEBI_23367, RDF.type, OWL.Class))

    # --- UniProt local UP_* classes ---
    for u in [n for n in nodes if isinstance(n, URIRef) and str(n).startswith(str(MAK)) and tail(n).startswith("UP_")]: # check new uniprot classes
        for _, _, o in list(outg.triples((u, RDFS.subClassOf, None))): # check what are their superclass
            outg.remove((u, RDFS.subClassOf, o)) # remove existing subclass relations
        outg.add((u, RDFS.subClassOf, PROTEIN)) # make them subclass of protein
        outg.add((u, RDF.type, OWL.Class)) # make them class

    # Protein ⊑ Target only
    for _, _, o in list(outg.triples((PROTEIN, RDFS.subClassOf, None))): # check superclasses of Protein class
        outg.remove((PROTEIN, RDFS.subClassOf, o)) # remove them
    outg.add((PROTEIN, RDFS.subClassOf, TARGET)) # make Protein class subclass of Target
    outg.add((PROTEIN, RDF.type, OWL.Class)) # Protein and Target are class
    outg.add((TARGET, RDF.type, OWL.Class))

    # --- Orphanet classes ---
    for ocls in [n for n in nodes if isinstance(n, URIRef) and str(n).startswith(str(ORDO)) and tail(n).startswith("Orphanet_")]: # we get all triples involving ORDO classes
        for _, _, o in list(outg.triples((ocls, RDFS.subClassOf, None))): # check their superclasses
            outg.remove((ocls, RDFS.subClassOf, o)) # remove them all
        if ocls != ORPHANET_C001: # if the current orpha class is not root orpha
            outg.add((ocls, RDFS.subClassOf, ORPHANET_C001)) # we make it subclass of root orpha
        outg.add((ocls, RDF.type, OWL.Class))
    for _, _, o in list(outg.triples((ORPHANET_C001, RDFS.subClassOf, None))):  # check superclasses of rooot orpha
        outg.remove((ORPHANET_C001, RDFS.subClassOf, o)) # we remove them all
    outg.add((ORPHANET_C001, RDFS.subClassOf, AUTOIMMUNE_DZ)) # we make root_orpha subclass of autoimmuneDisease
    outg.add((ORPHANET_C001, RDF.type, OWL.Class))
    outg.add((AUTOIMMUNE_DZ, RDF.type, OWL.Class))

    # --- UMLS CUI classes: classify by relations
    # Priority: disease over target; default Target if unseen.
    for c in cui_all:
        for _, _, o in list(outg.triples((c, RDFS.subClassOf, None))): # we check superclasses of all CUI classes
            outg.remove((c, RDFS.subClassOf, o))
        if c in cui_disease:
            outg.add((c, RDFS.subClassOf, AUTOIMMUNE_DZ)) # if it is invovled in the disease set, we make that CUI class a suclass of autoimmune disease
        elif c in cui_target:
            outg.add((c, RDFS.subClassOf, TARGET)) # if it is involved in the Target set, we make it subclass of Target
        else:
            outg.add((c, RDFS.subClassOf, TARGET))
        #outg.add((c, RDF.type, OWL.Class)) # if the CUI is not linked to disease or target, e don't add it anywhere
        pass

    # --- Remove Document and Relation classes entirely ---
    for banned in EXCLUDE_CLASSES:
        for t in list(outg.triples((banned, None, None))): # remove all triples where a banned class is subject or object, just to make sure
            outg.remove(t)
        for t in list(outg.triples((None, None, banned))):
            outg.remove(t)

force_constraints(out) # call function tha tapplies constraints to new KG

# Ensure CUI classes have labels (copy from source class first, then from *_instance)
def ensure_cui_labels(outg: Graph, gsrc: Graph):
    nodes = all_nodes(outg) # get all subject-objects from new KG
    cui_classes = [n for n in nodes
                   if isinstance(n, URIRef)
                   and "/concept/C" in str(n)
                   and not str(n).endswith("_instance")] # all CUI classes
    for c in cui_classes: # go through these classes
        if has_label(outg, c): # if it has label, skip
            continue
        for o in gsrc.objects(c, RDFS.label): # if that class has label in input graph, we add it
            outg.add((c, RDFS.label, o))
        for o in gsrc.objects(c, SKOS.prefLabel): # same for pre flabel
            outg.add((c, SKOS.prefLabel, o))
        if has_label(outg, c): # if has label now, skip
            continue
        inst = URIRef(str(c) + "_instance") # if still no label, we go look for labels in the instances of that class in original KG
        for o in gsrc.objects(inst, RDFS.label):
            outg.add((c, RDFS.label, o))
        for o in gsrc.objects(inst, SKOS.prefLabel):
            outg.add((c, SKOS.prefLabel, o))

ensure_cui_labels(out, g) # call the function tha tmakes sure CUI classes have labels


# 4) Serialize
out.serialize(destination=OUT_PATH, format="application/rdf+xml") # export to rdf file

print("Input triples:", len(g)) # display how many tripls we get
print("Output triples:", len(out))
print("Wrote:", OUT_PATH)
