import nibabel as nib
import numpy as np
import pandas as pd
import os
import warnings
import numpy as np
from nilearn import image
from itertools import combinations
from nilearn.maskers import NiftiMasker

thresholds = [round(0.5 + i * 0.05, 2) for i in range(11)]
thresholds_str = [f"{num:.2f}" for num in thresholds]
thresholds_percent = [str(int(num*100)) for num in thresholds]
cohorts = ['200_perm'+str(i) for i in range(1,6)]+['500_perm'+str(i) for i in range(1,6)]+['1000_perm'+str(i) for i in range(1,6)]+['3000_perm'+str(i) for i in range(1,6)]+['5000_perm'+str(i) for i in range(1,6)]+['7000_perm'+str(i) for i in range(1,6)]
dataPath = 'group/'
tracts = ['Cingulum_ALL','externalCapsule','ppn_thalamus']

# revise from https://github.com/demidenm/PyReliMRI/blob/main/pyrelimri/similarity.py
def image_similarity(imgfile1: str, imgfile2: str,
                     mask: str = None, thresh: float = None,
                     similarity_type: str = 'dice') -> float:
    # load list of images
    imagefiles = [imgfile1, imgfile2]
    img = [image.load_img(i) for i in imagefiles]
    assert img[0].shape == img[1].shape, 'images of different shape, ''image 1 {} and image 2 {}'.format(img[0].shape, img[1].shape)
    # mask image and binarize the image
    masker = NiftiMasker(mask_img=mask)
    imgdata = masker.fit_transform(img)
    imgdata[0, :] = (imgdata[0, :]>0).astype(int)
    imgdata[1, :] = (imgdata[1, :]>0).astype(int)
    if similarity_type.casefold() in ['dice']:
        # Intersection of images
        intersect = np.logical_and(imgdata[0, :], imgdata[1, :])
        dice_coeff = (2*intersect.sum()) / (float(imgdata[0, :].sum()+imgdata[1, :].sum()) + np.finfo(float).eps)
        coef = dice_coeff
    elif similarity_type.casefold() in ['jaccard']:
        intersect = np.logical_and(imgdata[0, :], imgdata[1, :])
        union = np.logical_or(imgdata[0, :], imgdata[1, :])
        jaccard_coeff = (intersect.sum()) / (float(union.sum()) + np.finfo(float).eps)
        coef = jaccard_coeff
    return coef


def pairwise_similarity(nii_filelist: list, mask: str = None,
                        thresh: float = None, similarity_type: str = 'Dice') -> pd.DataFrame:
    var_pairs = list(combinations(nii_filelist, 2))
    coef_df = pd.DataFrame(columns=['similar_coef', 'image_labels'])
    for img_comb in var_pairs:
        # select basename of file name(s)
        path = [os.path.basename(i) for i in img_comb]
        # calculate simiarlity
        val = image_similarity(imgfile1=img_comb[0], imgfile2=img_comb[1], mask=mask,
                               thresh=thresh, similarity_type=similarity_type)
        # for each pairwise come, save value + label to pandas df
        similarity_data = pd.DataFrame(np.column_stack((val, " ~ ".join([path[0], path[1]]))),
                                    columns=['similar_coef', 'image_labels'])
        coef_df = pd.concat([coef_df, similarity_data], axis=0, ignore_index=True)
    return coef_df


def calc_Dice(list_path, similarity_type):
	result = pairwise_similarity(nii_filelist=list_path,similarity_type=similarity_type)
	return result

result_cingulum = pd.DataFrame(columns=['similar_coef', 'image_labels'])
result_external = pd.DataFrame(columns=['similar_coef', 'image_labels'])
result_ppnthalamus = pd.DataFrame(columns=['similar_coef', 'image_labels'])
for i in range(len(thresholds)-1):
	data_cingulum = []
	data_external = []
	data_ppnthalamus = []
	for j in range(len(cohorts)):
		cohort = cohorts[j]
		print(f"{i} threshold {j} cohort.")
		Path = dataPath+'templates_b1000_ukb_threshold_'+thresholds_str[i]+'_'+cohort+'/NF/'
		cingulum = Path+'probtrack_'+tracts[0]+'_bilateral_consensus'+thresholds_percent[i]+'_template_ukb_threshold_'+thresholds_str[i]+'_'+cohort+'.nii.gz'
		external = Path+'probtrack_'+tracts[1]+'_bilateral_consensus'+thresholds_percent[i]+'_template_ukb_threshold_'+thresholds_str[i]+'_'+cohort+'.nii.gz'
		ppnthalamus = Path+'probtrack_'+tracts[2]+'_bilateral_consensus'+thresholds_percent[i]+'_template_ukb_threshold_'+thresholds_str[i]+'_'+cohort+'.nii.gz'
		data_cingulum.append(cingulum)
		data_external.append(external)
		data_ppnthalamus.append(ppnthalamus)
	result_cingulum = pd.concat([result_cingulum, calc_Dice(data_cingulum[0:len(cohorts)],'dice')],axis=0)
	result_external = pd.concat([result_external, calc_Dice(data_external[0:len(cohorts)],'dice')],axis=0)
	result_ppnthalamus = pd.concat([result_ppnthalamus, calc_Dice(data_ppnthalamus[0:len(cohorts)],'dice')],axis=0)

result_cingulum['image_labels'] = result_cingulum['image_labels'].str.replace('probtrack_Cingulum_ALL_bilateral_consensus', '')
result_cingulum['image_labels'] = result_cingulum['image_labels'].str.replace('_template_ukb_threshold', '')
result_external['image_labels'] = result_external['image_labels'].str.replace('probtrack_externalCapsule_bilateral_consensus', '')
result_external['image_labels'] = result_external['image_labels'].str.replace('_template_ukb_threshold', '')
result_ppnthalamus['image_labels'] = result_ppnthalamus['image_labels'].str.replace('probtrack_ppn_thalamus_bilateral_consensus', '')
result_ppnthalamus['image_labels'] = result_ppnthalamus['image_labels'].str.replace('_template_ukb_threshold', '')
result_cingulum.to_csv('result/template_robustness/cingulum_diceCoef.csv',index=False)
result_external.to_csv('result/template_robustness/external_diceCoef.csv',index=False)
result_ppnthalamus.to_csv('result/template_robustness/ppnthalamus_diceCoef.csv',index=False)