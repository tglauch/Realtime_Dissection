import time
import requests
from bs4 import BeautifulSoup
import datetime
import os
from requests.packages.urllib3.exceptions import InsecureRequestWarning
import smtplib

from email.mime.text import MIMEText
from astropy.time import Time
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

def check_for_new_gcn():
    print(datetime.datetime.now())
    with open('./gcns.txt', 'r') as f:
        gcns = f.read()
    all_gcn = gcns.split('\n')
    url = 'https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html'
    try:
        html= requests.get(url, verify=False,  timeout=10).text.encode('utf8')
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
            send_mail()
            ex_str = 'screen -d -m -S newevent bash -c \'python roi_analysis.py --gcn {}\''.format(full_href)
            print(ex_str)
            os.system(ex_str)
    with open('./gcns.txt', 'w+') as f:
        f.write('\n'.join(all_gcn))
    return

def send_mail():

    t = Time.now()
    ev_str = 'IC{}{:02d}{:02d}'.format(str(t.datetime.year)[-2:],
                                       t.datetime.month, t.datetime.day)
    ev_link='https://icecube.wisc.edu/~tglauch/Reports/{}/'.format(ev_str)

    # Create the container (outer) email message.
    msg = MIMEText('A new IceCube alert just came in...results will soon appear on {}'.format(ev_link))
    msg['Subject'] = 'New Alert {}'.format(ev_str)
    me = 'theo.glauch@tum.de'
    to = ['theo.glauch@tum.de', 'giommipaolo@gmail.com']
    msg['From'] = me
    msg['To'] = ', '.join(to)

    # Send the email via our own SMTP server.
    s = smtplib.SMTP('localhost')
    s.sendmail(me, to, msg.as_string())
    s.quit()

while True:
    check_for_new_gcn()
    time.sleep(60)
