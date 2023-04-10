from model import ExtractBindingSites

bin_sites = ExtractBindingSites().extract_binding_sites('3bj1', 'A')
print(bin_sites)