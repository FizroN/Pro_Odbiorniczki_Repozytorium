text = 'PawelKowalcze', Nbits = 4; % Zmieniając ilość bitów zmieniam ilość paczek w których jest wysyłany tekst
[numbers, bitsnum, bitschar] = text2numbers( text, Nbits ), pause
numbers2text( numbers, Nbits), pause 