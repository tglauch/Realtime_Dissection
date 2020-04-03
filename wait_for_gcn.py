import time
import requests
from bs4 import BeautifulSoup
import datetime
import os

def check_for_new_gcn():
    print(datetime.datetime.now())
    with open('./gcns.txt', 'r') as f:
        gcns = f.read()
    all_gcn = gcns.split('\n')
    url = 'https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html'
    try:
        html= requests.get(url).text.encode('utf8')
    except requests.exceptions.RequestException as e:
        print('Not able to information stuff from the GCN Website')
        return
    table = html.split('<!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX>')[1].split('<!--')[0]
    parsed_html = BeautifulSoup(table, "html5lib")
    for a in parsed_html.find_all('a', href=True):
        full_href  =  'https://gcn.gsfc.nasa.gov/' + a['href']
        if not full_href in all_gcn:
            print('NEW GCN')
            print(full_href)
            all_gcn.append(full_href)
            all_gcn = [i for i in all_gcn if i != '']
            os.system('screen -d -m -S newevent bash -c \'python roi_analysis.py --gcn {}\''.format(full_href))
    with open('./gcns.txt', 'w+') as f:
        f.write('\n'.join(all_gcn))
    return

while True:
    check_for_new_gcn()
    time.sleep(60)
