alphabet = 'abcdefghijklmnopqrstuvwxyz '
key = ''
text = ''

with open('key.txt', 'r') as infile:
    key = infile.readlines()[0].lower()

with open('plaintext.txt', 'r') as infile:
    text = infile.readlines()

with open('ciphertext.txt', 'w') as outfile:
    for line in text:
        keyIndices = [alphabet.index(k.lower()) for k in line]
        outfile.write(''.join(key[keyIndex] for keyIndex in keyIndices))
