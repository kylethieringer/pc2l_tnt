import os
import glob
import sleap.info.write_tracking_h5
from sleap.io.dataset import Labels

def main():
    exptList = glob.glob('/run/user/1000/gvfs/smb-share:server=cup.pni.princeton.edu,share=murthy/Kyle/flies/pc2l_tnt/**/**/*.proofread.slp')
    print(len(exptList))

    for filepath in exptList:
        print(os.path.basename(filepath))
        output_path = filepath.split('.slp')[0]+'.tracking.h5'
        if os.path.exists(output_path):
            continue

        video_callback = Labels.make_video_callback([os.path.dirname(filepath)])
        labels = Labels.load_file(filepath, video_search=video_callback)

        sleap.info.write_tracking_h5.main(labels, output_path=output_path, all_frames=True)



if __name__=="__main__":
    main()