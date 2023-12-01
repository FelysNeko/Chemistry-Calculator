def res(data, success):
    if success:
        return {'result': data, 'msg': '', 'success': True}
    else:
        return {'result': '', 'msg': data, 'success': False}
    