"""Small shared option parser for the focused application."""
def pairs(text):
    if not isinstance(text,str):raise ValueError('Column mapping and f-double-prime must be text')
    result={}
    if not text.strip():return result
    for item in text.split(','):
        if '=' not in item:raise ValueError('Use comma-separated name=value entries')
        key,value=map(str.strip,item.split('=',1))
        if not key or not value or key in result:raise ValueError('Mapping entries must have unique nonempty names and values')
        result[key]=value
    return result


from pathlib import Path


def read_sequences(path):
    records,identifier,parts = [],None,[]
    skip_description = False
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            if identifier is not None:
                records.append((identifier,''.join(parts)))
            identifier,parts = line[1:].strip(),[]
            skip_description = identifier[:3].upper() in ['P1;','F1;']
        elif skip_description:
            skip_description = False
        elif line.strip():
            if identifier is None:
                raise ValueError('Sequence must have a FASTA or PIR header')
            parts.append(''.join(line.split()).rstrip('*').upper())
    if identifier is not None:
        records.append((identifier,''.join(parts)))
    if not records or any(not seq or set(seq)-set('ACDEFGHIKLMNPQRSTVWY') for _,seq in records):
        raise ValueError('Invalid or empty protein sequence; supported codes are the 20 standard amino acids')
    return records
