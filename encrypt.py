alphabet = 'abcdefghijklmnopqrstuvwxyz '
key = ''
text = ''

with open('key.txt', 'r') as infile:
    key = infile.readlines()[0].lower()

with open('plaintext.txt', 'r') as infile:
    text = infile.readlines()

with open('ciphertext.txt', 'w') as outfile:
    for line in text:
        keyIndices = []
        for k in line:
            try:
                keyIndices.append(alphabet.index(k.lower()))
            except ValueError:
                print('Not found in alphabet or in key:', k.lower())
        outfile.write(''.join(key[keyIndex] for keyIndex in keyIndices))
