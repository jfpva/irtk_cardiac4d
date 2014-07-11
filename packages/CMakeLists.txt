SUBDIRS(transformation segmentation registration registration2 applications) #rview 

IF (BUILD_CARDIAC)
  SUBDIRS(cardiac)
ENDIF (BUILD_CARDIAC)


