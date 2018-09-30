from slackclient import SlackClient


def print_to_slack(text, pic=None):
    token = "xoxp-2661557430-30887140263-446963154231-478d91d13194c7fb8180734d000b9005"
    client_id = "2661557430.445727547252"
    client_secret = "c7446a6d5d5aeeea3f635914af240672"

    args = {'channel': "mario_tennis",
            'text': text}

    if pic is not None:
        args['attachments'] = [{"title": "Region around the Event",
                                "image_url": pic}]

    sc = SlackClient(token, client_id=client_id, client_secret=client_secret)
    x = sc.api_call("chat.postMessage", *args)
    return x
