#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function, division
import os, bz2, json, time
from datetime import timedelta
import cPickle as pickle

from pyhip import HIP

if __name__ == '__main__':
    # == == == == == == == == Part 1: Load ACTIVE dataset == == == == == == == == #
    # First time it gets loaded from the JSON format and writes essential fields into a pickle binary file.
    # check if the binary exists
    if not os.path.exists('./data/active-dataset.p'):
        print('>>> Converting ACTIVE dataset from JSON format to pickle... might take a while!')
        active_videos = {}
        with bz2.BZ2File('./data/active-dataset.json.bz2') as f:
            dataset = json.loads(f.readline())
            for video in dataset:
                active_videos[video['YoutubeID']] = (video['numShare'], video['dailyViewcount'], video['watchTime'])
        pickle.dump(active_videos, open('./data/active-dataset.p', 'wb'))

    print('>>> Loading the ACTIVE dataset from pickle...')
    active_videos = pickle.load(open('./data/active-dataset.p', 'rb'))

    # == == == == == == == == Part 2: Fit model and forecast future volume == == == == == == == == #
    start_time = time.time()
    # select exemplary video
    test_vid = 'X0ZEt_GZfkA'

    # uncomment next 2 lines if select random video
    # import random
    # test_vid = random.choice(active_videos.keys())

    print('>>> Fitting and forecasting for video: {0}'.format(test_vid))
    daily_share, daily_view, daily_watch = active_videos[test_vid]

    # uncomment next 2 lines if fit daily watch time
    # # convert watch time to hour unit
    # daily_watch = [x / 60 for x in daily_watch]

    num_train = 90
    num_test = 30
    num_initialization = 25

    hip_model = HIP()
    hip_model.initial(daily_share, daily_view, num_train, num_test, num_initialization)
    hip_model.fit_with_bfgs()
    hip_model.print_parameters()
    print('>>> Total fitting time: {0}'.format(str(timedelta(seconds=time.time() - start_time)))[:-3])
    hip_model.plot_func('YouTubeID={0}'.format(test_vid))
