import argparse
import math

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Create animation from nbody output.")
    sim_group = parser.add_mutually_exclusive_group(required=True)
    sim_group.add_argument('--galaxy', action='store_true', help="Select galaxy animation.")
    sim_group.add_argument('--solar', action='store_true', help="Select solar animation.")

    file_group = parser.add_mutually_exclusive_group(required=True)
    file_group.add_argument('--mp4', action='store_true', help="Save as mp4 file.")
    file_group.add_argument('--gif', action='store_true', help="Save as gif file.")

    return parser.parse_args()


def save_animation(ani, mp4=False):
    file_name = 'nbody_animation'

    print(f'Saving animation to {file_name} ...')

    fps= 1000 / ani.event_source.interval
    metadata = dict(title='n-body simulation', comment='Made with stdpar')

    if mp4:
        writer = animation.FFMpegWriter(
            fps=fps,
            metadata=metadata,
        )
        file_name += '.mp4'
    else:
        # gif writer
        writer = animation.PillowWriter(
            fps=fps,
            metadata=metadata,
        )
        file_name += '.gif'

    ani.save(
        file_name,
        writer=writer,
        savefig_kwargs={'pad_inches': 0},
    )


def read_points(file_name='positions.bin'):
    print(f'Reading {file_name}...')
    # read properties of file
    file_info = np.memmap(file_name, np.uint32, 'r', shape=3)

    # extract properties
    sim_size, steps, data_size = file_info
    dtype = np.float32 if data_size == 4 else np.float64 if data_size == 8 else 1 / 0

    # load data
    data = np.memmap(file_name, dtype, 'r', shape=(steps, 2, sim_size), offset=12)
    print(f'Loaded {data.shape}')

    return data


def animate_galaxy():
    data = read_points()

    # set up background
    fig, ax = plt.subplots(figsize=(6, 6))#, frameon=False)
    fig.tight_layout()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_axis_off()
    size = 500
    ax.set_xlim([-size, size])
    ax.set_ylim([-size, size])

    # create each frame
    artists = []
    for step in data[::10, ...]:
        n = step.shape[-1]
        step_1, step_2 = step[:, :n // 2], step[:, n // 2:]

        artist1 = ax.scatter(*step_1, marker='o', animated=True, color='red', s=1)
        artist2 = ax.scatter(*step_2, marker='o', animated=True, color='blue', s=1)
        artists.append([artist1, artist2])

    print(f'There are {len(artists)} frames')

    # build animation
    ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=100, blit=True, repeat_delay=1000)
    print('Animation created!')

    return ani


def animate_solar_system():
    # keep a frame a day
    data = read_points()[::24, ...]

    data_x = data[:, 0, :]
    data_y = data[:, 1, :]
    data_norm = np.sqrt(data_x ** 2 + data_y ** 2)
    multiplier = data_norm ** (1 / 15) / data_norm
    data_x = multiplier * data_x
    data_y = multiplier * data_y

    # planet info
    planets = ['sun', 'mecury', 'venus', 'earth', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    colours = ['yellow', 'grey', 'orange', 'blue', 'red', 'orange', 'orange', 'blue', 'blue']
    distances = list(range(len(planets)))
    size = 10

    # set up background
    fig, ax = plt.subplots(figsize=(6, 6))
    fig.tight_layout()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    ax.set_axis_off()
    ax.set_xlim([-size, size])
    ax.set_ylim([-size, size])

    # create each frame
    artists = []
    for xs, ys in zip(data_x, data_y):
        artist = ax.scatter(xs, ys, marker='o', color=colours)
        artists.append([artist])

    print(f'There are {len(artists)} frames')

    # build animation
    ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=20, blit=True)#, repeat_delay=1000)

    return ani





def main(args):
    if args.galaxy:
        ani = animate_galaxy()
    else:
        ani = animate_solar_system()

    save_animation(ani, mp4=args.mp4)


if __name__ == '__main__':
    main(parse_args())
