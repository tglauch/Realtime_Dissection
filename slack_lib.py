from slackclient import SlackClient


def print_to_slack(text, pic=None):
    token = "xoxp-2661557430-30887140263-447517735830-ec0a8c60dbebc0c4a0c87b8ffdc5c4d9"
    client_id = "2661557430.445727547252"
    client_secret = "c7446a6d5d5aeeea3f635914af240672"

    args = {'channel': 'test_channel',
            'text': text}

    sc = SlackClient(token, client_id=client_id, client_secret=client_secret)
    x = sc.api_call("chat.postMessage", **args)
    if pic is not None:
        y = sc.api_call('files.upload', file=open(pic, 'rb'), channels=['test_channel'])
    return x
